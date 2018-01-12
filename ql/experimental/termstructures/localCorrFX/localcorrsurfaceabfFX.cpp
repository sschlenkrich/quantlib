/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C)  2017 Cord Harms

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql\experimental\termstructures\localCorrFX\localcorrsurfaceabfFX.hpp>
#include <ql\experimental\termstructures\Helper\ParticleMethodUtils.hpp>

namespace QuantLib {

    LocalCorrSurfaceABFFX::LocalCorrSurfaceABFFX(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal)
    : LocalCorrSurfaceABF(processes, processToCal){
		corr0_ = RealStochasticProcess::MatA(2);
		corr1_ = RealStochasticProcess::MatA(2);
		
		for (size_t k = 0; k<2; ++k) corr0_[k].resize(2);
		for (size_t k = 0; k<2; ++k) corr1_[k].resize(2);
		
		//correlation is determined by lambda in abf-class, therefore simple correlation matrices

		corr0_[0][0] = 1;
		corr0_[1][0] = 0;
		corr0_[0][1] = 0;
		corr0_[1][1] = 1;

		corr1_[0][0] = 1;
		corr1_[1][0] = 1;
		corr1_[0][1] = 1;
		corr1_[1][1] = 1;
    }

	QuantLib::Real LocalCorrSurfaceABFFX::localFStrike(Time t, const RealStochasticProcess::VecA& X0) {
		QL_REQUIRE(X0.size() == 2, "Local Correlation for FX only works for two dimensional FX model.");

		return ParticleMethodUtils::getCrossFX(processes_[0]->x0() * std::exp(X0[0]) , (processes_[1]->x0() * std::exp(X0[1])));
		
	}
	
	void LocalCorrSurfaceABFFX::accept(AcyclicVisitor& v) {
		Visitor<LocalCorrSurfaceABFFX>* v1 =
			dynamic_cast<Visitor<LocalCorrSurfaceABFFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABF::accept(v);
	}

	Real LocalCorrSurfaceABFFX::localCorrImplTeq0(Time t, const RealStochasticProcess::VecA& X0, bool extrapolate) {
		
		//smiled surface will through an error, therefore assume one minute ahead
		t = 1.0 / (365 * 24 * 60);
		
		Real s1 = processes_[0]->x0() * std::exp(X0[0]);
		Real s2 = processes_[1]->x0() * std::exp(X0[1]);
		Real vol1 = processes_[0]->localVolatility()->localVol(t, s1, extrapolate);
		Real vol2 = processes_[1]->localVolatility()->localVol(t, s2, extrapolate);
		Real vol3 = processToCal_->localVolatility()->localVol(t, ParticleMethodUtils::getCrossFX(s1 , s2), extrapolate);

		return (vol1*vol1 + vol2*vol2 - vol3*vol3) / (2 * vol1*vol2);
	}
	
}

