/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Sebastian Schlenkrich
*/

/*! \file vanillalocalvoltermstructures.cpp
    \brief swaption volatility term structure based on VanillaLocalVolModel
*/


#include <ql/termstructures/volatility/smilesection.hpp>
#include <ql/indexes/swapindex.hpp>

#include <ql/experimental/templatemodels/vanillalocalvol/vanillalocalvoltermstructures.hpp>
#include <ql/experimental/templatemodels/vanillalocalvol/vanillalocalvolsmilesection.hpp>


namespace QuantLib {

	VanillaLocalVolSwaptionVTS::VanillaLocalVolSwaptionVTS(
		const boost::shared_ptr<SwaptionVolatilityStructure>&                                     atmVolTS,
		const std::vector< std::vector< boost::shared_ptr<VanillaLocalVolModelSmileSection> > >&  smiles,
		const std::vector< Period >&                                                              swapTerms,
		const boost::shared_ptr<SwapIndex>&                                                       index)
	: SwaptionVolatilityStructure(atmVolTS->referenceDate(),atmVolTS->calendar(),atmVolTS->businessDayConvention(), atmVolTS->dayCounter(), atmVolTS->volatilityType()),
		atmVolTS_(atmVolTS), smiles_(smiles), swapTerms_(swapTerms), index_(index) {
		QL_REQUIRE(atmVolTS_, "atmVolTS required");
		QL_REQUIRE(smiles_.size() == swapTerms_.size(), "smiles_.size()==swapTerms_.size() required");
		for (size_t k = 1; k < swapTerms_.size(); ++k) {
			QL_REQUIRE(months(swapTerms_[k - 1]) < months(swapTerms_[k]), "months(swapTerms_[k-1])<months(swapTerms_[k]) required");
		}
		for (size_t k = 0; k < smiles_.size(); ++k) {
			QL_REQUIRE(smiles_[k].size() > 0, "smiles_[k].size()>0 required");
			for (size_t i = 1; i < smiles_[k].size(); ++i) {
				QL_REQUIRE(smiles_[k][i - 1]->exerciseTime() < smiles_[k][i]->exerciseTime(), "smiles_[k][i-1]->exerciseTime()<smiles_[k][i]->exerciseTime() required");
			}
		}
	}

	boost::shared_ptr<SmileSection> VanillaLocalVolSwaptionVTS::smileSectionImpl( Time optionTime, Time swapLength) const {
		if (smiles_.size() == 0) return atmVolTS_->smileSection(optionTime, swapLength);  // fall back
		// we need to convert back times to periods and dates
		Period optionDays((Integer)round(365.0 * optionTime), Days);
		Date optionDate = referenceDate() + optionDays;
		Period swapTerm((Integer)round(12.0 * swapLength), Months);
		// first we interpolate in expiry direction for two swap term columns
		std::vector< size_t > idxj(2);  // swap term indices enclosing swap length
		std::vector< boost::shared_ptr<VanillaLocalVolModelSmileSection> > smilesj(2);
		idxj[0] = swapTerms_.size() - 1;
		idxj[1] = 0;
		while ((idxj[0] > 0)                     && (months(swapTerms_[idxj[0]]) > months(swapTerm))) --idxj[0];
		while ((idxj[1] < swapTerms_.size() - 1) && (months(swapTerms_[idxj[1]]) < months(swapTerm))) ++idxj[1];
		if (smiles_[idxj[0]].size() == 0) return atmVolTS_->smileSection(optionTime, swapLength);  // fall back
		if (smiles_[idxj[1]].size() == 0) return atmVolTS_->smileSection(optionTime, swapLength);  // fall back
		Real rhoj = 0.5;
		if (months(swapTerms_[idxj[0]]) < months(swapTerms_[idxj[1]]))
			rhoj = (round(12.0 * swapLength) - months(swapTerms_[idxj[0]])) / (months(swapTerms_[idxj[1]]) - months(swapTerms_[idxj[0]]));
		for (size_t k = 0; k < 2; ++k) {
			std::vector< size_t > idxi(2);  // expiry indices enclosing optionTime
			idxi[0] = smiles_[idxj[k]].size() - 1;
			idxi[1] = 0;
			while ((idxi[0] > 0)                           && (smiles_[idxj[k]][idxi[0]]->exerciseTime() > optionTime)) --idxi[0];
			while ((idxi[1] < smiles_[idxj[k]].size() - 1) && (smiles_[idxj[k]][idxi[1]]->exerciseTime() < optionTime)) ++idxi[1];
			Real rhoi = 0.5;
			if (smiles_[idxj[k]][idxi[0]]->exerciseTime() < smiles_[idxj[k]][idxi[1]]->exerciseTime())
				rhoi = (optionTime - smiles_[idxj[k]][idxi[0]]->exerciseTime()) / (smiles_[idxj[k]][idxi[1]]->exerciseTime() - smiles_[idxj[k]][idxi[0]]->exerciseTime());
			// now we may interpolate in expiry direction
			Rate forward = index_->clone(swapTerms_[idxj[k]])->fixing(optionDate);  // this might fail if optionDate is no good business day
			Real atmVol = atmVolTS_->volatility(optionTime, swapLength, forward);
			smilesj[k] = boost::shared_ptr<VanillaLocalVolModelSmileSection>(
				new VanillaLocalVolModelSmileSection(optionDate, forward, atmVol, smiles_[idxj[k]][idxi[0]], smiles_[idxj[k]][idxi[1]], rhoi, true, dayCounter(), referenceDate(), VolatilityType(), atmVolTS_->shift(optionDate, swapLength)));				
		}
		Rate forward = index_->clone(swapTerm)->fixing(optionDate);  // this might fail if optionDate is no good business day
		Real atmVol = atmVolTS_->volatility(optionTime, swapLength, forward); // we want to be as accuarate as possible on the ATM interpolation thus using times
		boost::shared_ptr<VanillaLocalVolModelSmileSection> smile(new VanillaLocalVolModelSmileSection(optionDate, forward, atmVol, smilesj[0], smilesj[1], rhoj, true, dayCounter(), referenceDate(), VolatilityType(), atmVolTS_->shift(optionDate, swapLength)));
		return smile;
	}

}

