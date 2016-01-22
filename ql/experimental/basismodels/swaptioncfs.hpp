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

	class IborLegCashFlows {
    protected:
		Date                                 refDate_;   // today, base for time calculations w.r.t. Act/365 (Fixed) 
        Leg                                  floatLeg_;
		std::vector<Real>                    floatTimes_;
		std::vector<Real>                    floatWeights_;
	public:
		inline const Leg&               floatLeg()     const { return floatLeg_;       }
		inline const std::vector<Real>& floatTimes()   const { return floatTimes_;     }
		inline const std::vector<Real>& floatWeights() const { return  floatWeights_;  }
        IborLegCashFlows  ( const Leg&                        iborLeg,
			                const Handle<YieldTermStructure>& discountCurve,
						    bool                              contTenorSpread = true );
	};


    class SwapCashFlows : public IborLegCashFlows {
    protected:
		// resulting cash flows as leg
        Leg                                  fixedLeg_;
		std::vector<Real>                    fixedTimes_;
		std::vector<Real>                    fixedWeights_;
		std::vector<Real>                    annuityWeights_;
	public:
		SwapCashFlows ( const boost::shared_ptr<VanillaSwap>& swap,
			            const Handle<YieldTermStructure>&     discountCurve,
						bool                                  contTenorSpread = true );
        // inspectors
        inline const Leg&               fixedLeg()         const  { return fixedLeg_; }
		inline const std::vector<Real>& fixedTimes()	   const  { return  fixedTimes_;	 }
		inline const std::vector<Real>& fixedWeights()     const  { return  fixedWeights_;	 }
		inline const std::vector<Real>& annuityWeights()   const  { return  annuityWeights_; }		
    };


    class SwaptionCashFlows : public SwapCashFlows {
	protected:
		boost::shared_ptr<Swaption>          swaption_;
		std::vector<Real>                    exerciseTimes_;
	public:
   	    SwaptionCashFlows ( const boost::shared_ptr<Swaption>&    swaption,
			                const Handle<YieldTermStructure>&     discountCurve,
						    bool                                  contTenorSpread = true );
        // inspectors
		inline const boost::shared_ptr<Swaption> swaption() const { return swaption_; }
		inline const std::vector<Real>& exerciseTimes()    const  { return  exerciseTimes_;  }
	};


}

#endif /* #ifndef quantlib_SwaptionCashFlows_hpp */
