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



namespace QuantLib {

	// Base class for template payoffs
    template <class DateType, class PassiveType, class ActiveType>
	class MCPayoffT {
	protected:
		typedef MCSimulationT<DateType, PassiveType, ActiveType>        SimulationType;
		typedef typename MCSimulationT<DateType, PassiveType, ActiveType>::Path  PathType;

	    DateType observationTime_;
	public:
	    MCPayoffT( const DateType observationTime ) : observationTime_(observationTime) { }
		// inspectors
		inline DateType observationTime() { return observationTime_; }
		// generic payoff(observationTime, p) needs to be implemented by derived classes
        inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) = 0;
		// discounted payoff for NPV valuation
		inline virtual ActiveType discountedAt(const boost::shared_ptr<PathType>& p) { return at(p) / p->numeraire(observationTime_); }

		// generic pricer
		class Pricer {
		protected:
			std::vector< boost::shared_ptr<MCPayoffT> > payoffs_;
			boost::shared_ptr<SimulationType>           simulation_;
		public:

			inline static std::vector<ActiveType> at( const boost::shared_ptr<MCPayoffT>&       payoff,
				                                      const boost::shared_ptr<SimulationType>&  simulation) {
			    std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k=0; k<simulation->nPaths(); ++k) res[k] = payoff->at(simulation->path(k));
				return res;
			}

			inline static std::vector<ActiveType> discountedAt( const boost::shared_ptr<MCPayoffT>&       payoff,
				                                                const boost::shared_ptr<SimulationType>&  simulation) {
			    std::vector<ActiveType> res(simulation->nPaths());
				for (size_t k=0; k<simulation->nPaths(); ++k) res[k] = payoff->discountedAt(simulation->path(k));
				return res;
			}

			inline static ActiveType NPV(const std::vector< boost::shared_ptr<MCPayoffT> >&  payoffs, 
				                         const boost::shared_ptr<SimulationType>&            simulation) {
				ActiveType npv = 0.0;
				for (size_t k=0; k<payoffs.size(); ++k) {
					ActiveType npvk = 0;
					for (size_t n=0; n<simulation->nPaths(); ++n) {
						npvk += payoffs[k]->discountedAt(simulation->path(n));
					}
					npv += npvk;
				}
				return npv / simulation->nPaths();
			}

			Pricer (const std::vector< boost::shared_ptr<MCPayoffT> >&   payoffs, 
				    const boost::shared_ptr<SimulationType>&             simulation)
					: payoffs_(payoffs), simulation_(simulation) { }

			inline ActiveType NPV() {  return NPV(payoffs_,simulation_);  }
		};


		// basic payoffs and operations

		// a deterministic flow known in advance (undiscounted)
        class FixedAmount : public MCPayoffT {
		protected:
			ActiveType amount_;
		public:
			FixedAmount( const ActiveType amount ) : MCPayoffT(0.0), amount_(amount) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return amount_; }
		};

		// a x + y  (undiscounted)
        class Axpy : public MCPayoffT {
		protected:
			ActiveType a_;
			boost::shared_ptr<MCPayoffT> x_, y_;
		public:
			Axpy( const ActiveType                      a,
				  const boost::shared_ptr<MCPayoffT>&   x,
				  const boost::shared_ptr<MCPayoffT>&   y )
				  : MCPayoffT(0.0), a_(a), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				ActiveType res = a_ * x_->at(p);
				if (y_) res += y_->at(p);
				return res;
			}
		};

		// x * y  (undiscounted)		
        class Mult : public MCPayoffT {
		
		protected:
			boost::shared_ptr<MCPayoffT> x_, y_;
        
		public:
			Mult ( const boost::shared_ptr<MCPayoffT>&   x,
				   const boost::shared_ptr<MCPayoffT>&   y )
				  : MCPayoffT(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return x_->at(p) * y_->at(p);
			}
		};

		// x / y  (undiscounted)		
		class Division : public MCPayoffT {

		protected:
			boost::shared_ptr<MCPayoffT> x_, y_;

		public:
			Division(const boost::shared_ptr<MCPayoffT>&   x,
				     const boost::shared_ptr<MCPayoffT>&   y)
				    : MCPayoffT(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return x_->at(p) / y_->at(p);
			}
		};

		// max{x,y}  (undiscounted)
        class Max : public MCPayoffT {
		protected:
			boost::shared_ptr<MCPayoffT> x_, y_;
		public:
			Max ( const boost::shared_ptr<MCPayoffT>&   x,
				  const boost::shared_ptr<MCPayoffT>&   y )
				  : MCPayoffT(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return (x_->at(p)>y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
		};

		// min{x,y}  (undiscounted)
        class Min : public MCPayoffT {
		protected:
			boost::shared_ptr<MCPayoffT> x_, y_;
		public:
			Min ( const boost::shared_ptr<MCPayoffT>&   x,
				  const boost::shared_ptr<MCPayoffT>&   y )
				  : MCPayoffT(0.0), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return (x_->at(p)<y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
		};

		// (re-)set paydate of a payoff (for discounting)
        class Pay : public MCPayoffT {
		protected:
			boost::shared_ptr<MCPayoffT> x_;
		public:
			Pay ( const boost::shared_ptr<MCPayoffT>&   x,
				  const DateType                        payTime )
				  : MCPayoffT(payTime), x_(x) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return x_->at(p);
			}
		};

		// simple discounted cash payment
		class Cash : public MCPayoffT {
		protected:
			DateType payTime_;
		public:
			Cash( DateType obsTime, DateType payTime ) : MCPayoffT(obsTime), payTime_(payTime) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				//if (payTime_<=observationTime()) return (ActiveType)1.0;
				return p->zeroBond(observationTime(),payTime_);  // catch any exception in path, simulation or model
			}
		};

		// 1 unit of modelled asset
		class Asset : public MCPayoffT {
		protected:
			std::string alias_;
		public:
			Asset(DateType obsTime, const std::string alias) : MCPayoffT(obsTime), alias_(alias) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return p->asset(observationTime(), alias_); }
		};

		// 1 unit call or put exercised at observation time and settled at pay time
		class VanillaOption : public MCPayoffT {
		protected:
			DateType    payTime_;
			PassiveType callOrPut_;
			PassiveType strike_;
		public:
			VanillaOption( DateType obsTime, DateType payTime, PassiveType strike, PassiveType callOrPut ) : MCPayoffT(obsTime), payTime_(payTime), strike_(strike), callOrPut_(callOrPut) { }
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (payTime_<observationTime()) return (ActiveType)0.0;
				ActiveType DF = p->zeroBond(observationTime(),payTime_);
				ActiveType S  = p->asset(observationTime(),"");
				ActiveType V  = callOrPut_ * DF * (S - strike_);
				return (V>0.0) ? (V) : ((ActiveType)0.0);
			}
		};

		// cache result in case it is requested repeatedly
		class Cache : public MCPayoffT {
		protected:
			boost::shared_ptr<MCPayoffT> x_;
			boost::shared_ptr<PathType> lastPath_;
			ActiveType                  lastPayoff_;
		public:
			Cache(const boost::shared_ptr<MCPayoffT>&   x) : MCPayoffT(x->observationTime()), x_(x) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (lastPath_!=p) {
					lastPath_ = p;
					lastPayoff_ = x_->at(p);
				}
				return lastPayoff_;
			}
		};

		// logical operators
		class Logical : public MCPayoffT {
		protected:
			// seems we need to wrap oparators to construct a pointer to it
			static bool equal       (const ActiveType& x, const ActiveType& y) { return (x == y); }
			static bool notEqual    (const ActiveType& x, const ActiveType& y) { return (x != y); }
			static bool less        (const ActiveType& x, const ActiveType& y) { return (x < y);  }
			static bool lessEqual   (const ActiveType& x, const ActiveType& y) { return (x <= y); }
			static bool greater     (const ActiveType& x, const ActiveType& y) { return (x > y);  }
			static bool greaterEqual(const ActiveType& x, const ActiveType& y) { return (x >= y); }
			static bool and         (const ActiveType& x, const ActiveType& y) { return ((x!=(ActiveType)0.0)&&(y!=(ActiveType)0.0)); }
			static bool or          (const ActiveType& x, const ActiveType& y) { return ((x!=(ActiveType)0.0)||(y!=(ActiveType)0.0)); }
			// this is the actual pointer to the operator
			bool(*op_)(const ActiveType&, const ActiveType&);
			// these are the operands
			boost::shared_ptr<MCPayoffT> x_, y_;
		public:
			Logical(const boost::shared_ptr<MCPayoffT>&   x,
				    const boost::shared_ptr<MCPayoffT>&   y,
				    const std::string&                    op) : MCPayoffT(0.0), x_(x), y_(y) {
				op_ = &equal; // this is a very bad default
				if (op == "==") op_ = &equal;
				if (op == "!=") op_ = &notEqual;
				if (op == "<")  op_ = &less;
				if (op == "<=") op_ = &lessEqual;
				if (op == ">")  op_ = &greater;
				if (op == ">=") op_ = &greaterEqual;
				if (op == "&&") op_ = &and;
				if (op == "||") op_ = &or;
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if ((*op_)(x_->at(p), y_->at(p))) return (ActiveType)(1.0);
				else return (ActiveType)(0.0);
			}
		};

		class IfThenElse : public MCPayoffT {
		protected:
			boost::shared_ptr<MCPayoffT> x_, y_, z_;
		public:
			IfThenElse(const boost::shared_ptr<MCPayoffT>&   x,
				const boost::shared_ptr<MCPayoffT>&   y,
				const boost::shared_ptr<MCPayoffT>&   z) : MCPayoffT(0.0), x_(x), y_(y), z_(z) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				if (x_->at(p) != ((ActiveType)0.0)) return y_->at(p);
				return z_->at(p);
			}
		};

	};



}

#endif  /* ifndef quantlib_templatemcpayoff_hpp */
