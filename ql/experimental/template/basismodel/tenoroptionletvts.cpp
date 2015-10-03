/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Sebastian Schlenkrich
*/

/*! \file tenoroptionletvts.cpp
    \brief caplet volatility term structure based on volatility transformation
*/

#include <ql/experimental/template/basismodel/tenoroptionletvts.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/dategenerationrule.hpp>
#include <ql/math/rounding.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>
#include <ql/exercise.hpp>


namespace QuantLib {

	TenorOptionletVTS::TenorOptionletVTS(
		                   const Handle<OptionletVolatilityStructure>&            baseVTS,
			               const boost::shared_ptr<IborIndex>&                    baseIndex,
		                   const boost::shared_ptr<IborIndex>&                    targIndex,
						   const boost::shared_ptr<CorrelationStructure>&         correlation)
						   : OptionletVolatilityStructure(baseVTS->referenceDate(),baseVTS->calendar(),baseVTS->businessDayConvention(),baseVTS->dayCounter()),
		                  baseVTS_(baseVTS), baseIndex_(baseIndex), targIndex_(targIndex), correlation_(correlation) {
	    // check that target frequency is a multiple of base frequency
	}


	TenorOptionletVTS::TenorOptionletSmileSection::TenorOptionletSmileSection(
		    const TenorOptionletVTS& volTS,
			const Time               optionTime )
			: SmileSection(optionTime,volTS.baseVTS_->dayCounter(),Normal,0.0), correlation_(volTS.correlation_) {
		// we assume that long (target) tenor is a multiple of short (base) tenor
		// first we need the long tenor start and end date
		Real oneDayAsYear  = volTS.dayCounter().yearFraction(volTS.referenceDate(),volTS.referenceDate()+1);
		Date exerciseDate  = volTS.referenceDate() + ((BigInteger)ClosestRounding(0)(optionTime/oneDayAsYear));
		Date effectiveDate = volTS.baseIndex_->fixingCalendar().advance(exerciseDate, volTS.baseIndex_->fixingDays()*Days);
		Date maturityDate  = volTS.baseIndex_->fixingCalendar().advance(effectiveDate,volTS.targIndex_->tenor(),Unadjusted,false);
		// now we can set up the short tenor schedule
		Schedule baseFloatSchedule(effectiveDate, maturityDate, volTS.baseIndex_->tenor(),volTS.baseIndex_->fixingCalendar(),ModifiedFollowing,Unadjusted,DateGeneration::Backward,false);
        // set up scalar attributes
		fraRateTarg_ = volTS.targIndex_->fixing(exerciseDate);
		Time yfTarg   = volTS.targIndex_->dayCounter().yearFraction(effectiveDate,maturityDate);
		for (Size k=0; k<baseFloatSchedule.dates().size()-1; ++k) {
			Date startDate  = baseFloatSchedule.dates()[k];
			Date fixingDate = volTS.baseIndex_->fixingCalendar().advance(startDate,(-1*volTS.baseIndex_->fixingDays())*Days);
			Time yearFrac   = volTS.baseIndex_->dayCounter().yearFraction(baseFloatSchedule.dates()[k],baseFloatSchedule.dates()[k+1]);
			// set up vector attributes
			baseSmileSection_.push_back(volTS.baseVTS_->smileSection(fixingDate,true));
			startTimeBase_.push_back(volTS.dayCounter().yearFraction(volTS.referenceDate(),startDate));
			fraRateBase_.push_back(volTS.baseIndex_->fixing(fixingDate));
			v_.push_back( yearFrac / yfTarg * (1.0 + yfTarg*fraRateTarg_) / (1.0 + yearFrac*fraRateBase_[k]) );
		}	    
	}

    Volatility TenorOptionletVTS::TenorOptionletSmileSection::volatilityImpl(Rate strike) const {
		Real sum_v = 0.0;
		for (Size k=0; k<v_.size(); ++k) sum_v += v_[k];
		std::vector<Real> volBase(v_.size());
		for (Size k=0; k<fraRateBase_.size(); ++k) {
			Real strike_k = (strike - (fraRateTarg_ - sum_v*fraRateBase_[k])) / sum_v;
			volBase[k] = baseSmileSection_[k]->volatility(strike_k,Normal,0.0);
		}
		Real var=0.0;
		for (Size i=0; i<volBase.size(); ++i) {
			var += v_[i]*v_[i]*volBase[i]*volBase[i];
			for (Size j=i+1; j<volBase.size(); ++j) {
				Real corr = (*correlation_)(startTimeBase_[i],startTimeBase_[j]);
				var += 2.0*corr*v_[i]*v_[j]*volBase[i]*volBase[j];
			}
		}
		Real vol = sqrt(var);
		return vol;
	}


}

