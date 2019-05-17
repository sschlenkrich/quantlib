/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcpayoff.hpp
    \brief generic payoff interface for MC simulation
	
*/


#ifndef quantlib_templatemcpayoff_hpp
#define quantlib_templatemcpayoff_hpp


#include <ql/experimental/templatemodels/montecarlo/mcsimulationT.hpp>
#include <algorithm>
#include <vector>



namespace QuantLib {

	// Base class for template payoffs
    template <class DateType, class PassiveType, class ActiveType>
	class MCPayoffT {
	protected:
		typedef MCSimulationT<DateType, PassiveType, ActiveType>        SimulationType;
		typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

		DateType observationTime_;
	public:
		MCPayoffT(const DateType observationTime) : observationTime_(observationTime) { }
		// inspectors
		inline DateType observationTime() { return observationTime_; }
		// calculate observation times recursively
		inline virtual std::set<DateType> observationTimes() { std::set<DateType> s; s.insert(observationTime_); return s; }
		// generic payoff(observationTime, p) needs to be implemented by derived classes
		virtual ActiveType at(const boost::shared_ptr<PathType>& p) = 0;
		// discounted payoff for NPV valuation
		inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p) / p->numeraire(observationTime_); }
		// return a clone but with changed observation time; this effectively allows considering a payoff as an index
		inline virtual boost::shared_ptr<MCPayoffT> at(const DateType t) { QL_FAIL("at(t) not implemented."); }

		// for convenience, derive a union of two sets and return a new set
		static std::set<DateType> unionTimes(const std::set<DateType> s1, const std::set<DateType> s2) {
			std::set<DateType> s(s1);
			s.insert(s2.begin(), s2.end());
			return s;
		}

		// generic pricer
		class Pricer {
		protected:
			std::vector< boost::shared_ptr<MCPayoffT> > payoffs_;
			boost::shared_ptr<SimulationType>           simulation_;
		public:

			inline static std::vector<ActiveType> at(const boost::shared_ptr<MCPayoffT>&       payoff,
				const boost::shared_ptr<SimulationType>&  simulation) {
				std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k = 0; k < simulation->nPaths(); ++k) res[k] = payoff->at(simulation->path(k));
				return res;
			}

			inline static std::vector<ActiveType> discountedAt(const boost::shared_ptr<MCPayoffT>&       payoff,
				const boost::shared_ptr<SimulationType>&  simulation) {
				std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k = 0; k < simulation->nPaths(); ++k) res[k] = payoff->discountedAt(simulation->path(k));
				return res;
			}

			inline static ActiveType NPV(const std::vector< boost::shared_ptr<MCPayoffT> >&  payoffs,
				const boost::shared_ptr<SimulationType>&            simulation) {
				ActiveType npv = 0.0;
				for (size_t k = 0; k < payoffs.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n = 0; n < simulation->nPaths(); ++n) {
						npvk += payoffs[k]->discountedAt(simulation->path(n));
					}
					npv += npvk;
				}
				return npv / simulation->nPaths();
			}

			inline static std::vector<ActiveType> NPVs(const std::vector< boost::shared_ptr<MCPayoffT> >&  payoffs,
				                                       const boost::shared_ptr<SimulationType>&            simulation) {
				std::vector<ActiveType> res(payoffs.size());
				for (size_t n = 0; n < simulation->nPaths(); ++n) {
					for (size_t k = 0; k < payoffs.size(); ++k) {
						res[k] += payoffs[k]->discountedAt(simulation->path(n));
					}
				}
				for (size_t k = 0; k < payoffs.size(); ++k) {
					res[k] /= simulation->nPaths();
				}
				return res;
			}

			Pricer(const std::vector< boost::shared_ptr<MCPayoffT> >&   payoffs,
				const boost::shared_ptr<SimulationType>&             simulation)
				: payoffs_(payoffs), simulation_(simulation) { }

			inline ActiveType NPV() { return NPV(payoffs_, simulation_); }
		};
	};

	// Base template payoffs
    template <class DateType, class PassiveType, class ActiveType>
	class BasePayoffT {
		// basic payoffs and operations
	protected:
	    typedef MCPayoffT<DateType, PassiveType, ActiveType>  PayoffType;
		typedef MCSimulationT<DateType, PassiveType, ActiveType>  SimulationType;
		typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

	public:

		class Clone : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_;
		public:
			Clone(const boost::shared_ptr<PayoffType>&   x,
				  const DateType                        observationTime) : PayoffType(observationTime), x_(x->at(observationTime)) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return x_->at(p); }
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Clone(x_, t)); }
			inline virtual std::set<DateType> observationTimes() { return x_->observationTimes(); }
		};

		// a deterministic flow known in advance (undiscounted)
		class FixedAmount : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			ActiveType amount_;
		public:
			FixedAmount(const ActiveType amount) : PayoffType(0.0), amount_(amount) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return amount_; }
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new FixedAmount(amount_)); }
		};

		// (re-)set paydate of a payoff (for discounting)
		class Pay : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_;
		public:
			Pay(const boost::shared_ptr<PayoffType>&   x,
				const DateType                        payTime)
				: PayoffType(payTime), x_(x) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return x_->at(p); }
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Pay(x_->at(t), PayoffType::observationTime())); }
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(PayoffType::observationTimes(), x_->observationTimes()); }
		};

		// simple discounted cash payment
		class Cash : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			DateType payTime_;
		public:
			Cash(DateType obsTime, DateType payTime) : PayoffType(obsTime), payTime_(payTime) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				//if (payTime_<=observationTime()) return (ActiveType)1.0;
				return p->zeroBond(PayoffType::observationTime(), payTime_);  // catch any exception in path, simulation or model
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Cash(t, payTime_)); }

		};

		// 1 unit of modelled asset
		class Asset : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			std::string alias_;
		public:
			Asset(DateType obsTime, const std::string alias) : PayoffType(obsTime), alias_(alias) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return p->asset(PayoffType::observationTime(), alias_); }
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Asset(t, alias_)); }
		};

		// return the continuous barrier no-hit probability
		class AssetBarrierNoHit : public MCPayoffT<DateType,PassiveType,ActiveType> {
			std::string alias_;
			DateType tStart_, tEnd_;
			PassiveType downBarrier_, upBarrier_;
			PassiveType downOrUpOrBoth_; // down (-1), up (+1), both (0)
		public:
			AssetBarrierNoHit(DateType tStart, DateType tEnd, PassiveType downBarrier, PassiveType upBarrier, PassiveType downOrUpOrBoth, const std::string alias)
				: PayoffType(tEnd), tStart_(tStart), tEnd_(tEnd), alias_(alias),
				downBarrier_(downBarrier), upBarrier_(upBarrier), downOrUpOrBoth_(downOrUpOrBoth) {
				QL_REQUIRE(tStart_ < tEnd_, "AssetBarrierNoHit: tStart < tEnd required.");
				QL_REQUIRE(downBarrier_ < upBarrier_, "AssetBarrierNoHit: downBarrier < upBarrier required.");
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return p->assetBarrierNoHit(tStart_, tEnd_, downBarrier_, upBarrier_, downOrUpOrBoth_, alias_); }
			inline virtual std::set<DateType> observationTimes() { std::set<DateType> s; s.insert(tStart_); s.insert(tEnd_); return s; }
		};


		// 1 unit call or put exercised and settled at observation time
		class VanillaOption : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			std::string alias_;
			PassiveType callOrPut_;
			PassiveType strike_;
		public:
			VanillaOption(DateType obsTime, const std::string alias, PassiveType strike, PassiveType callOrPut) : PayoffType(obsTime), alias_(alias), strike_(strike), callOrPut_(callOrPut) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType S = p->asset(PayoffType::observationTime(), alias_);
				ActiveType V = callOrPut_ * (S - strike_);
				return (V > 0.0) ? (V) : ((ActiveType)0.0);
			}
		};

		// cache result in case it is requested repeatedly
		class Cache : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_;
			boost::shared_ptr<PathType> lastPath_;
			ActiveType                  lastPayoff_;
		public:
			Cache(const boost::shared_ptr<PayoffType>&   x) : PayoffType(x->observationTime()), x_(x) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (lastPath_ != p) {
					lastPath_ = p;
					lastPayoff_ = x_->at(p);
				}
				return lastPayoff_;
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Cache(x_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { return x_->observationTimes(); }
		};

		// arithmetics and functions applied to payoffs

		// a x + y  (undiscounted)
		class Axpy : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			ActiveType a_;
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Axpy(const ActiveType                      a,
				const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y)
				: PayoffType(0.0), a_(a), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				ActiveType res = a_ * x_->at(p);
				if (y_) res += y_->at(p);
				return res;
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Axpy(a_, x_->at(t), y_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { 
				return (y_)?(PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes())):(x_->observationTimes());
			}
		};

		// x * y  (undiscounted)		
		class Mult : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Mult(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y)
				: PayoffType(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return x_->at(p) * y_->at(p);
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Mult(x_->at(t), y_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()); }
		};

		// x / y  (undiscounted)		
		class Division : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Division(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y)
				: PayoffType(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return x_->at(p) / y_->at(p);
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Division(x_->at(t), y_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()); }
		};

		// max{x,y}  (undiscounted)
		class Max : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Max(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y)
				: PayoffType(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return (x_->at(p) > y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Max(x_->at(t), y_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()); }
		};

		// min{x,y}  (undiscounted)
		class Min : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Min(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y)
				: PayoffType(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return (x_->at(p) < y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) { return boost::shared_ptr<PayoffType>(new Min(x_->at(t), y_->at(t))); }
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()); }
		};

		// logical operators
		class Logical : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			// seems we need to wrap oparators to construct a pointer to it
			static bool equal(const ActiveType& x, const ActiveType& y) { return (x == y); }
			static bool notEqual(const ActiveType& x, const ActiveType& y) { return (x != y); }
			static bool less(const ActiveType& x, const ActiveType& y) { return (x < y); }
			static bool lessEqual(const ActiveType& x, const ActiveType& y) { return (x <= y); }
			static bool greater(const ActiveType& x, const ActiveType& y) { return (x > y); }
			static bool greaterEqual(const ActiveType& x, const ActiveType& y) { return (x >= y); }
			static bool and_(const ActiveType& x, const ActiveType& y) { return ((x != (ActiveType)0.0) && (y != (ActiveType)0.0)); }
			static bool or_(const ActiveType& x, const ActiveType& y)  { return ((x != (ActiveType)0.0) || (y != (ActiveType)0.0)); }
			// this is the actual pointer to the operator
			bool(*op_)(const ActiveType&, const ActiveType&);
			// these are the operands
			boost::shared_ptr<PayoffType> x_, y_;
		public:
			Logical(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y,
				const std::string&                    op) : PayoffType(0.0), x_(x), y_(y) {
				op_ = &equal; // this is a very bad default
				if (op == "==") op_ = &equal;
				if (op == "!=") op_ = &notEqual;
				if (op == "<")  op_ = &less;
				if (op == "<=") op_ = &lessEqual;
				if (op == ">")  op_ = &greater;
				if (op == ">=") op_ = &greaterEqual;
				if (op == "&&") op_ = &and_;
				if (op == "||") op_ = &or_;
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if ((*op_)(x_->at(p), y_->at(p))) return (ActiveType)(1.0);
				else return (ActiveType)(0.0);
			}
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()); }
		};

		class IfThenElse : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			boost::shared_ptr<PayoffType> x_, y_, z_;
		public:
			IfThenElse(const boost::shared_ptr<PayoffType>&   x,
				const boost::shared_ptr<PayoffType>&   y,
				const boost::shared_ptr<PayoffType>&   z) : PayoffType(0.0), x_(x), y_(y), z_(z) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (x_->at(p) != ((ActiveType)0.0)) return y_->at(p);
				return z_->at(p);
			}
			inline virtual std::set<DateType> observationTimes() { return PayoffType::unionTimes(PayoffType::unionTimes(x_->observationTimes(), y_->observationTimes()), z_->observationTimes()); }
		};

		class Basket : public MCPayoffT<DateType,PassiveType,ActiveType> {
		protected:
			std::vector< boost::shared_ptr<PayoffType> > underlyings_;
			std::vector< PassiveType > weights_;
			bool rainbow_;
			class Descending {
				const boost::shared_ptr<PathType>& p_;
			public:
				Descending(const boost::shared_ptr<PathType>& p) : p_(p) {}
				inline bool operator() (const boost::shared_ptr<PayoffType>& x, const boost::shared_ptr<PayoffType>& y) {
					return (x->at(p_) > y->at(p_));
				}
			};
		public:
			Basket(const std::vector<boost::shared_ptr<PayoffType>>& underlyings,
				const std::vector< PassiveType >                 weights,
				bool                                             rainbow)
				: PayoffType(0.0), underlyings_(underlyings), weights_(weights), rainbow_(rainbow) {
				QL_REQUIRE(underlyings_.size() > 0, "Basket underlyings required");
				QL_REQUIRE(underlyings_.size() == weights_.size(), "Basket dimension mismatch");
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (rainbow_) std::sort(underlyings_.begin(), underlyings_.end(), Descending(p));
				ActiveType res = 0;
				for (size_t k = 0; k < underlyings_.size(); ++k) res += weights_[k] * underlyings_[k]->at(p);
				return res;
			}
			inline virtual boost::shared_ptr<PayoffType> at(const DateType t) {
				std::vector<boost::shared_ptr<PayoffType>> underlyingsAt;
				for (size_t k = 0; k < underlyings_.size(); ++k) underlyingsAt.push_back(underlyings_[k]->at(t));
				return boost::shared_ptr<PayoffType>(new Basket(underlyingsAt, weights_, rainbow_));
			}
			inline virtual std::set<DateType> observationTimes() {
				std::set<DateType> s;
				for (size_t k = 0; k < underlyings_.size(); ++k) s = PayoffType::unionTimes(s, underlyings_[k]->observationTimes());
				return s;
			}
		};

	};
}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
