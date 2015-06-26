/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file swaptioncfs.hpp
    \brief translate swaption into deterministic fixed and float cash flows
*/

#ifndef quantlib_swaptioncfs_hpp
#define quantlib_swaptioncfs_hpp

#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/option.hpp>

namespace QuantLib {

    class SwaptionCashFlows {
    protected:
		Date                                 refDate_;   // today, base for time calculations w.r.t. Act/365 (Fixed) 
		boost::shared_ptr<Swaption>          swaption_;  // reference to underlying swaption
		// resulting cash flows as leg
        Leg                                  fixedLeg_;
        Leg                                  floatLeg_;
		// resultng cash flows as raw data
		std::vector<Real>                    exerciseTimes_;
		std::vector<Real>                    fixedTimes_;
		std::vector<Real>                    floatTimes_;
		std::vector<Real>                    fixedWeights_;
		std::vector<Real>                    floatWeights_;
	public:
		// constructor to map a swaption to bond option according to spread model
		SwaptionCashFlows ( const boost::shared_ptr<Swaption>& swaption,
			                const Handle<YieldTermStructure>& discountCurve,
						    bool                              contTenorSpread = true );
        // inspectors
		inline const boost::shared_ptr<Swaption> swaption() const { return swaption_; }
        inline const Leg&                        fixedLeg() const { return fixedLeg_; }
		inline const Leg&                        floatLeg() const { return floatLeg_; }
		// assemble cash flow values and pay times w.r.t. yield curve
		inline const std::vector<Real>& exerciseTimes()   { return  exerciseTimes_; }
		inline const std::vector<Real>& fixedTimes()	  { return 	fixedTimes_;	}
		inline const std::vector<Real>& floatTimes()	  { return 	floatTimes_;	}
		inline const std::vector<Real>& fixedWeights()    { return 	fixedWeights_;	}
		inline const std::vector<Real>& floatWeights()    { return 	floatWeights_;	}
		
    };

}

#endif /* #ifndef quantlib_SwaptionCashFlows_hpp */
