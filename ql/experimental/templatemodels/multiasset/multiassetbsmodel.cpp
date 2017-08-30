/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Sebastian Schlenkrich

*/


#include <ql/experimental/templatemodels/auxilliaries/svdT.hpp>

#include <ql/experimental/templatemodels/multiasset/multiassetbsmodel.hpp>


namespace QuantLib {

	MultiAssetBSModel::MultiAssetBSModel(
		const Handle<YieldTermStructure>&                                               termStructure,
		const std::vector<std::string>&                                                 aliases,
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const RealStochasticProcess::MatA&                                              correlations)
	: termStructure_(termStructure), processes_(processes) {
		QL_REQUIRE(processes_.size() > 0, "No BS processes supplied");
		QL_REQUIRE(processes_.size() == aliases.size(), "Number of processes doesn't match aliases");
		for (size_t k = 0; k < aliases.size(); ++k) index_[aliases[k]] = k;
		DT_ = TemplateAuxilliaries::svdSqrt(correlations);
	}

	// initial values for simulation
	inline RealStochasticProcess::VecP MultiAssetBSModel::initialValues() {
		return RealStochasticProcess::VecP(size(), 0.0);
	}
	// a[t,X(t)]
	inline RealStochasticProcess::VecA MultiAssetBSModel::drift(const QuantLib::Time t, const VecA& X) {
		RealStochasticProcess::VecA nu(processes_.size());
		// todo: make sure all processes use same domestic/riskfree rate...
		for (Size k = 0; k < processes_.size(); ++k) {
			QuantLib::Real S = processes_[k]->x0() * std::exp(X[k]);
			nu[k] = processes_[k]->drift(t, S);
		}
		return nu;
	}
	// b[t,X(t)]
	inline RealStochasticProcess::MatA MultiAssetBSModel::diffusion(const QuantLib::Time t, const VecA& X) {
		RealStochasticProcess::MatA b(DT_);
		for (Size i = 0; i < size(); ++i) {
			QuantLib::Real S = processes_[i]->x0() * std::exp(X[i]);
			QuantLib::Real sigma = processes_[i]->diffusion(t, S);
			for (Size j = 0; j < factors(); ++j) {
				b[i][j] *= sigma;
			}
		}
		return b;
	}

	inline void MultiAssetBSModel::evolve(const QuantLib::Time t0, const VecA& X0, const QuantLib::Time dt, const VecD& dW, VecA& X1) {
		// we approximate the local vol diffusion by a Black-Scholes process on [t, t+dt]
		for (Size i = 0; i < X1.size(); ++i) {
			X1[i] = 0.0;
			for (Size j = 0; j < dW.size(); ++j) X1[i] += DT_[i][j] * dW[j];
			// sigma represents the average volatility on [t, t+dt]
			// here we use a first very simple approximation
			Real S = processes_[i]->x0() * std::exp(X0[i]);
			Real sigma = processes_[i]->diffusion(t0, S);
			// We may integrate the drift exactly (given the approximate volatility)
			Real B_d = processes_[i]->riskFreeRate()->discount(t0) / processes_[i]->riskFreeRate()->discount(t0 + dt);
			Real B_f = processes_[i]->dividendYield()->discount(t0) / processes_[i]->dividendYield()->discount(t0 + dt);
			X1[i] = X0[i] + std::log(B_d / B_f) - 0.5 * sigma * sigma * dt + sigma * X1[i] * std::sqrt(dt);
		}
	}

}



