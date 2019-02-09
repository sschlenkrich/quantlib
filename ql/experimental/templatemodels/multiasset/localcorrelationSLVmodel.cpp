/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019, Cord Harms

*/


#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>

#include <ql/experimental/templatemodels/multiasset/localCorrelationslvmodel.hpp>


namespace QuantLib {

	LocalCorrelationSLVModel::LocalCorrelationSLVModel(
		const Handle<YieldTermStructure>&                                               termStructure,
		const std::vector<std::string>&                                                 aliases,
		const std::vector<boost::shared_ptr<QuantLib::HestonSLVProcess>>& processes,
		const Handle<LocalCorrTermStructure>& localCorrTermStructure)
	: MultiAssetSLVModel(termStructure,aliases,processes), 
		localCorrTermStructure_(localCorrTermStructure) {
		corrMatrix_ = RealStochasticProcess::MatA(processes.size());

		for (size_t k = 0; k<processes.size(); ++k) corrMatrix_[k].resize(processes.size());
	}

	inline void LocalCorrelationSLVModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
		// we merely need to calculate the correlation matrix
		localCorrTermStructure_->localCorr(corrMatrix_, t0, X0 , true); //true because first X0 are always 0 in simulation.
		//DT_ = TemplateAuxilliaries::svdSqrt(corrMatrix_);
		TemplateAuxilliaries::performCholesky(corrMatrix_, corrMatrix_.size());
		DT_ = corrMatrix_;
		//now we call the function of multiassetmodel
		MultiAssetSLVModel::evolve(t0, X0, dt, dW, X1);
	}

}



