/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatemcpayoff.hpp
    \brief generic payoff interface for MC simulation
	
*/


#ifndef quantlib_templatemcswap_hpp
#define quantlib_templatemcswap_hpp


#include <ql/experimental/template/montecarlo/templatemcpayoff.hpp>
#include <ql/experimental/template/basismodel/swaptioncfs.hpp>

#include <ql/indexes/iborindex.hpp>
#include <ql/indexes/swapindex.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>


namespace QuantLib {

	template <class DateType, class PassiveType, class ActiveType>
	class TemplateMC {
	protected:
		typedef TemplateMCSimulation<DateType, PassiveType, ActiveType>        SimulationType;
		typedef typename TemplateMCSimulation<DateType, PassiveType, ActiveType>::Path  PathType;

	public:
	// basic payoffs based on Libor and swap rates

        class FixedAmount : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			ActiveType amount_;
		public:
			FixedAmount( const ActiveType amount ) : TemplateMCPayoff(0.0), amount_(amount) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return amount_; }
		};

        class LiborRate : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			DateType    fixingTime_, startTime_, endTime_;
			PassiveType oneOverDaycount_;
			PassiveType D_; // tenor basis
		public:
			LiborRate( const DateType                        fixingTime,
				       const DateType                        startTime,
					   const DateType                        endTime,
					   const boost::shared_ptr<IborIndex>&   iborIndex,
					   const Handle<YieldTermStructure>&     discYTS )
					   : TemplateMCPayoff(fixingTime), fixingTime_(fixingTime), startTime_(startTime), endTime_(endTime) {
				Date today      = discYTS->referenceDate(); // check if this is the correct date...
				Date fixingDate = today + ((BigInteger)ClosestRounding(0)(fixingTime_*365.0)); // assuming act/365 day counting
				Date startDate  = today + ((BigInteger)ClosestRounding(0)(startTime_ *365.0)); // assuming act/365 day counting
				Date endDate    = today + ((BigInteger)ClosestRounding(0)(endTime_   *365.0)); // assuming act/365 day counting
				PassiveType liborForward = iborIndex->fixing(fixingDate,true);
				PassiveType daycount     = iborIndex->dayCounter().yearFraction(startDate,endDate);
				oneOverDaycount_ = 1.0/daycount;
				// tenor basis calculation
				D_ = (1.0 + daycount*liborForward) * discYTS->discount(endDate) / discYTS->discount(startDate);
			}			        
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) {
				return ( p->zeroBond(fixingTime_,startTime_) / p->zeroBond(fixingTime_,endTime_) * D_ - 1.0 ) * oneOverDaycount_;
			}
	    };

		class SwapRate : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			boost::shared_ptr<TemplateMCPayoff<DateType, PassiveType, ActiveType>::SwapRate > swaprate_;
		public:
			SwapRate( const DateType                        fixingTime,
					  const boost::shared_ptr<SwapIndex>&   swapIndex,
					  const Handle<YieldTermStructure>&     discYTS    ) : TemplateMCPayoff(fixingTime) {
				Date today      = discYTS->referenceDate(); // check if this is the correct date...
				Date fixingDate = today + ((BigInteger)ClosestRounding(0)(fixingTime*365.0)); // assuming act/365 day counting
				SwapCashFlows scf(swapIndex->underlyingSwap(fixingDate),discYTS,true);  // assume continuous tenor spreads
				swaprate_ = boost::shared_ptr<TemplateMCPayoff::SwapRate>(new TemplateMCPayoff::SwapRate(
					fixingTime, scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights() ) );
			}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { return swaprate_->at(p); }
		};

        class Axpy : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			ActiveType a_;
			boost::shared_ptr<TemplateMCPayoff> x_, y_;
		public:
			Axpy( const ActiveType                             a,
				  const boost::shared_ptr<TemplateMCPayoff>&   x,
				  const boost::shared_ptr<TemplateMCPayoff>&   y )
				  : TemplateMCPayoff(x->observationTime()), a_(a), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				ActiveType res = a_ * x_->at(p);
				if (y_) res += y_->at(p);
				return res;
			}
		};

        class Max : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			boost::shared_ptr<TemplateMCPayoff> x_, y_;
		public:
			Max ( const boost::shared_ptr<TemplateMCPayoff>&   x,
				  const boost::shared_ptr<TemplateMCPayoff>&   y )
				  : TemplateMCPayoff(x->observationTime()), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return (x_->at(p)>y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
		};

        class Min : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			boost::shared_ptr<TemplateMCPayoff> x_, y_;
		public:
			Min ( const boost::shared_ptr<TemplateMCPayoff>&   x,
				  const boost::shared_ptr<TemplateMCPayoff>&   y )
				  : TemplateMCPayoff(x->observationTime()), x_(x), y_(y) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return (x_->at(p)<y_->at(p)) ? (x_->at(p)) : (y_->at(p));
			}
		};

        class Pay : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			boost::shared_ptr<TemplateMCPayoff> x_;
		public:
			Pay ( const boost::shared_ptr<TemplateMCPayoff>&   x,
				  const DateType                               payTime )
				  : TemplateMCPayoff(payTime), x_(x) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return x_->at(p);
			}
		};


	    // coupon decorating a payoff with start and pay date for organisation in legs
        class Coupon : public TemplateMCPayoff<DateType, PassiveType, ActiveType> {
		protected:
			boost::shared_ptr<TemplateMCPayoff> x_;
			DateType                            startTime_, payTime_;
		public:
			Coupon ( const boost::shared_ptr<TemplateMCPayoff>&   x,
				     const DateType                               startTime,
				     const DateType                               payTime )
		        : TemplateMCPayoff(payTime), x_(x), startTime_(startTime), payTime_(payTime) {}
			Coupon ( const boost::shared_ptr<TemplateMCPayoff>&   x )
				: TemplateMCPayoff(x->observationTime()), x_(x), startTime_(x->observationTime()), payTime_(x->observationTime()) {}
			inline virtual ActiveType at(const boost::shared_ptr<PathType>& p) { 
				return x_->at(p);
			}
			// inspectors
			inline const DateType startTime() const { return startTime_; }
			inline const DateType payTime()   const { return payTime_;   }
		};

	    // a coupon leg as a ordered list of coupons
		class Leg : public std::vector< boost::shared_ptr<Coupon> > {
		public:
			static bool firstIsLess( const boost::shared_ptr<Coupon>& first, const boost::shared_ptr<Coupon>& second ) { return first->startTime() < second->startTime(); } 
			Leg ( const std::vector< boost::shared_ptr<Coupon> >&   coupons ) : std::vector< boost::shared_ptr<Coupon> >(coupons.begin(),coupons.end())  {
				// check that all coupons are consistent
				// sort by start time
				std::sort(this->begin(),this->end(),firstIsLess );
			}
		};

	    // a swap as a set of coupon legs (e.g. structured, funding, notional exchanges)
		class Swap : public std::vector< boost::shared_ptr<Leg> > {
		public:
			Swap ( const std::vector< boost::shared_ptr<Leg> >&   legs ) : std::vector< boost::shared_ptr<Leg> >(legs.begin(),legs.end())  { }
		};


		// this is the key structure for AMC valuation
		class CancellableNote {
		private:
			// underlying
			std::vector< boost::shared_ptr<Leg> >  underlyings_;       // the underlying coupon legs
			// call features
			std::vector< DateType >                callTimes_;            // exercise times
			std::vector< boost::shared_ptr<Leg> >  earlyRedemptions_;     // strikes payed at exercise
			std::vector< boost::shared_ptr<Leg> >  regressionVariables_;  // regression variables at exercise
		public:
			CancellableNote ( const std::vector< boost::shared_ptr<Leg> >&  underlyings,
				              const std::vector< DateType >&                callTimes,
							  const std::vector< boost::shared_ptr<Leg> >&  earlyRedemptions,
							  const std::vector< boost::shared_ptr<Leg> >&  regressionVariables )
							  : underlyings_(underlyings), callTimes_(callTimes), earlyRedemptions_(earlyRedemptions), regressionVariables_(regressionVariables) {
			    // sanity checks
			}
			// inspectors
			inline const std::vector< boost::shared_ptr<Leg> >& underlyings()         const { return underlyings_;         }
			inline const std::vector< DateType >&               callTimes()           const { return callTimes_;           }
			inline const std::vector< boost::shared_ptr<Leg> >& earlyRedemptions()    const { return earlyRedemptions_;    }
			inline const std::vector< boost::shared_ptr<Leg> >& regressionVariables() const { return regressionVariables_; }
		};


	}; // class TemplateMC

}

#endif  /* ifndef quantlib_templatemcswap_hpp */ 
