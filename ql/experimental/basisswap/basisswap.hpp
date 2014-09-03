/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswap.hpp
    \brief (Cross currency) Interest rate swap consisting of several legs
*/

#ifndef quantlib_basisswap_hpp
#define quantlib_basisswap_hpp

#include <ql/instrument.hpp>
#include <ql/cashflow.hpp>
#include <ql/instruments/swap.hpp>



namespace QuantLib {

	// class BasisSwap extends Swap pricing by FX rate for foreign currency legs
	// and par rate valuation for botstrapping via RateHelper
    class BasisSwap : public Swap {
    protected:
        // data members for par rate/spread valuation
		Size parLegIndex_;     // define the leg number on which par rate should be evaluated
		bool calcParSpread_;   // define if par spread (true) or par rate (false) should be evaluated
    public:
		class arguments;
		class results;
		class engine;
        //! \name Constructors
        /*! Multi leg constructor. */
        BasisSwap(const std::vector<Leg>&  legs,
                  const std::vector<bool>& payer,
                  Size                     parLegIndex   = 0,
				  bool                     calcParSpread = 0
				  ) : Swap(legs, payer), parLegIndex_(parLegIndex), calcParSpread_(calcParSpread) {
			  QL_REQUIRE(parLegIndex_<legs_.size(), "Wrong par leg index");
		}
		Real fairRate() const;
		// Inspectors
		const std::vector<Leg>&  legs()  const { return legs_;          }
		const std::vector<bool> payer() const {
			std::vector<bool> boolpayer(payer_.size());
			for (Size i=0; i<boolpayer.size(); ++i) boolpayer[i] = (payer_[i]<0) ? true : false;
			return boolpayer;
		}
		Size  parLegIndex()              const { return parLegIndex_;   }
		bool  calcParSpread()            const { return calcParSpread_; }
    };

	class BasisSwap::arguments : public Swap::arguments { };

    class BasisSwap::results : public Swap::results { };

    class BasisSwap::engine : public GenericEngine<BasisSwap::arguments,
                                                   BasisSwap::results> {};


}

#endif
