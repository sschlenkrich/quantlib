/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Sebastian Schlenkrich

*/



#ifndef quantlib_multiassetbsmodel_hpp
#define quantlib_multiassetbsmodel_hpp


#include <ql/experimental/templatemodels/stochasticprocessT.hpp>

#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {


	// We model a multi-asset local volatility model by means of the normalized log-processes X_i = log[S_i/S_(0)]

	class MultiAssetBSModel : public RealStochasticProcess {
	private:
		Handle<YieldTermStructure>                                               termStructure_;  // domestic discounting term structure
		std::map<std::string, size_t>                                            index_;
		std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>> processes_;
		RealStochasticProcess::MatA                                              DT_;  // D^T D = Correlations
	public:
		MultiAssetBSModel(const Handle<YieldTermStructure>&                                               termStructure,
			              const std::vector<std::string>&                                                 aliases,
			              const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
			              const RealStochasticProcess::MatA&                                              correlations);

		// dimension of X
		inline virtual size_t size() { return processes_.size(); }
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline virtual size_t factors() { return processes_.size(); }
		// initial values for simulation
		inline virtual RealStochasticProcess::VecP initialValues();
		// a[t,X(t)]
		inline virtual RealStochasticProcess::VecA drift(const QuantLib::Time t, const VecA& X);
		// b[t,X(t)]
		inline virtual RealStochasticProcess::MatA diffusion(const QuantLib::Time t, const VecA& X);

		inline virtual void evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1);

		// default implementation, zero interest rates
		inline virtual QuantLib::Real numeraire(const QuantLib::Time t, const VecA& X) { return 1.0/processes_[0]->riskFreeRate()->discount(t); }

		// default implementation, zero interest rates
		inline virtual QuantLib::Real zeroBond(const QuantLib::Time t, const QuantLib::Time T, const VecA& X) { return processes_[0]->riskFreeRate()->discount(T) / processes_[0]->riskFreeRate()->discount(t); }

		// default implementation for single-asset models
		inline virtual QuantLib::Real asset(const QuantLib::Time t, const VecA& X) { return processes_[0]->x0() * std::exp(X[0]); }

		// default implementation for single-asset models
		inline virtual QuantLib::Real asset(const QuantLib::Time t, const VecA& X, const std::string& alias) { return processes_[index_.at(alias)]->x0() * std::exp(X[index_.at(alias)]); }

	};

}



#endif  /* ifndef quantlib_multiassetbsmodel_hpp */
