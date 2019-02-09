/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Cord Harms

*/



#ifndef quantlib_localcorrelationbsmodel_hpp
#define quantlib_localcorrelationbsmodel_hpp


#include <ql/experimental/templatemodels/multiasset/multiassetbsmodel.hpp>
#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {


	// We model a multi-asset local volatility model by means of the normalized log-processes X_i = log[S_i/S_(0)]

	class LocalCorrelationBSModel : public MultiAssetBSModel {
	private:
		RealStochasticProcess::MatA corrMatrix_;
		Handle<LocalCorrTermStructure> localCorrTermStructure_;
	public:
		LocalCorrelationBSModel(const Handle<YieldTermStructure>&                                         termStructure,
			              const std::vector<std::string>&                                                 aliases,
			  			  const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
						  const Handle<LocalCorrTermStructure>&											  localCorrTermStructure);

		inline virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);
	};

}



#endif  /* ifndef quantlib_localcorrelationbsmodel_hpp */
