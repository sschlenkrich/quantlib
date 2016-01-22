/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich
*/

/*! \file tenorswaptionvts.hpp
    \brief swaption volatility term structure based on volatility transformation
*/

#ifndef quantlib_tenorswaptionvts_hpp
#define quantlib_tenorswaptionvts_hpp

#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/option.hpp>

namespace QuantLib {

    class TenorSwaptionVTS : public SwaptionVolatilityStructure {
    protected:

		class TenorSwaptionSmileSection : public SmileSection {
        protected:
            boost::shared_ptr<SmileSection> baseSmileSection_;
			Real    swapRateBase_;
			Real    swapRateTarg_;
			Real    swapRateFinl_;
			Real    lambda_;
			Real    annuityScaling_;
			// implement transformation formula
            virtual Volatility volatilityImpl(Rate strike) const;
		public:
			// constructor includes actual transformation details
			TenorSwaptionSmileSection ( const TenorSwaptionVTS& volTS,
				                        Time                    optionTime,
                                        Time                    swapLength);

		    // further SmileSection interface methods
			virtual Real minStrike() const { return baseSmileSection_->minStrike()+swapRateTarg_-swapRateBase_; }
			virtual Real maxStrike() const { return baseSmileSection_->maxStrike()+swapRateTarg_-swapRateBase_; }
			virtual Real atmLevel() const  { return swapRateFinl_;  }
		};

		Handle<SwaptionVolatilityStructure>            baseVTS_;
		Handle<YieldTermStructure>                     discountCurve_;

		boost::shared_ptr<IborIndex>                   baseIndex_;
		boost::shared_ptr<IborIndex>                   targIndex_;
		Period                                         baseFixedFreq_;
		Period                                         targFixedFreq_;
		DayCounter                                     baseFixedDC_;
		DayCounter                                     targFixedDC_;

	public:
		// constructor
		TenorSwaptionVTS( const Handle<SwaptionVolatilityStructure>&            baseVTS,
		                  const Handle<YieldTermStructure>&                     discountCurve,
		                  const boost::shared_ptr<IborIndex>&                   baseIndex,
		                  const boost::shared_ptr<IborIndex>&                   targIndex,
		                  const Period&                                         baseFixedFreq,
		                  const Period&                                         targFixedFreq,
		                  const DayCounter&                                     baseFixedDC,
		                  const DayCounter&                                     targFixedDC)
						  : SwaptionVolatilityStructure(baseVTS->referenceDate(),baseVTS->calendar(),baseVTS->businessDayConvention(),baseVTS->dayCounter()),
		                  baseVTS_(baseVTS), discountCurve_(discountCurve), baseIndex_(baseIndex), targIndex_(targIndex),
						  baseFixedDC_(baseFixedDC), targFixedDC_(targFixedDC) { }

		// Termstructure interface

		//! the latest date for which the curve can return values
		virtual Date maxDate() const { return baseVTS_->maxDate(); }

		// SwaptionVolatility interface

		//! the minimum strike for which the term structure can return vols
		virtual Rate minStrike() const { return baseVTS_->minStrike(); }
        //! the maximum strike for which the term structure can return vols
        virtual Rate maxStrike() const {return baseVTS_->maxStrike(); }


		// SwaptionVolatilityStructure interface

		//! the largest length for which the term structure can return vols
		virtual const Period& maxSwapTenor() const { return baseVTS_->maxSwapTenor(); }

        virtual boost::shared_ptr<SmileSection> smileSectionImpl(
                                                Time optionTime,
												Time swapLength) const { return boost::shared_ptr<SmileSection>(new TenorSwaptionSmileSection(*this,optionTime,swapLength) ); }

        virtual Volatility volatilityImpl(Time optionTime,
                                          Time swapLength,
										  Rate strike) const { return smileSectionImpl(optionTime,swapLength)->volatility(strike,VolatilityType::Normal,0.0); }


    };

}

#endif /* #ifndef quantlib_tenorswaptionvts_hpp */
