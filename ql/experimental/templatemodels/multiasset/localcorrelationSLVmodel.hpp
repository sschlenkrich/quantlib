/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2019, Cord Harms

*/



#ifndef quantlib_localcorrelationslvmodel_hpp
#define quantlib_localcorrelationslvmodel_hpp


#include <ql/experimental/templatemodels/multiasset/multiassetslvmodel.hpp>
#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/experimental/processes/hestonslvprocess.hpp>

namespace QuantLib {


	// We model a multi-asset local volatility model by means of the normalized processes S_i

	class LocalCorrelationSLVModel : public MultiAssetSLVModel {
	private:
		RealStochasticProcess::MatA corrMatrix_;
		Handle<LocalCorrTermStructure> localCorrTermStructure_;
	public:
		LocalCorrelationSLVModel(const Handle<YieldTermStructure>&                                         termStructure,
			              const std::vector<std::string>&                                                 aliases,
			  			  const std::vector<boost::shared_ptr<QuantLib::HestonSLVProcess>>&				  processes,
						  const Handle<LocalCorrTermStructure>&											  localCorrTermStructure);

		inline virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);
	};

}



#endif  /* ifndef quantlib_localcorrelationslvmodel_hpp */
