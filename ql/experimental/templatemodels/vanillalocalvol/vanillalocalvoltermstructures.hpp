/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich
*/

/*! \file vanillalocalvoltermstructures.hpp
    \brief swaption volatility term structure based on VanillaLocalVolModel
*/

#ifndef quantlib_vanillalocalvoltermstructures_hpp
#define quantlib_vanillalocalvoltermstructures_hpp


#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>


namespace QuantLib {

	class VanillaLocalVolModelSmileSection;
	class SwapIndex;

    class VanillaLocalVolSwaptionVTS : public SwaptionVolatilityStructure {
    private:
		const boost::shared_ptr<SwaptionVolatilityStructure> atmVolTS_;
		const std::vector< std::vector< boost::shared_ptr<VanillaLocalVolModelSmileSection> > > smiles_;
		const std::vector< Period > swapTerms_;
		const boost::shared_ptr<SwapIndex> index_;  // a template for all the swap index per swap terms
	protected:
		virtual boost::shared_ptr<SmileSection> smileSectionImpl(
			Time optionTime,
			Time swapLength) const;

		virtual Volatility volatilityImpl(Time optionTime,
			Time swapLength,
			Rate strike) const {
			return smileSectionImpl(optionTime, swapLength)->volatility(strike);
		}

    public:
		VanillaLocalVolSwaptionVTS(
			const boost::shared_ptr<SwaptionVolatilityStructure>&                                     atmVolTS,
			const std::vector< std::vector< boost::shared_ptr<VanillaLocalVolModelSmileSection> > >&  smiles,
			const std::vector< Period >&                                                              swapTerms,
			const boost::shared_ptr<SwapIndex>&                                                       index);

        virtual const Date&   referenceDate() const { return atmVolTS_->referenceDate(); }
        virtual const Period& maxSwapTenor()  const { return atmVolTS_->maxSwapTenor();  }
		virtual       Date    maxDate()       const { return atmVolTS_->maxDate();       }
		virtual       Rate    minStrike()     const { return atmVolTS_->minStrike();     }
		virtual       Rate    maxStrike()     const { return atmVolTS_->maxStrike();     }

    };


}

#endif
