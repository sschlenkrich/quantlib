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

namespace QuantLib {

    LocalCorrSurfaceABFFX::LocalCorrSurfaceABFFX(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&  		    processToCal,
		boost::shared_ptr<CalibratorLocalCorrInt>&										calibratorLocalCorr)
    : LocalCorrSurfaceABF(processes, processToCal), calibratorLocalCorr_(calibratorLocalCorr){
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

	QuantLib::Real LocalCorrSurfaceABFFX::localF(Time t, const RealStochasticProcess::VecA& X0,
		bool extrapolate) const {
		QL_REQUIRE(X0.size()==2,"Local Correlation for FX only works for two dimensional FX model.")
		if (X0[1] !=0 ) return interpolatorF_(t, X0[0]/X0[1],extrapolate);
		return interpolatorF_(t, 0, extrapolate);
	}
	
	void LocalCorrSurfaceABFFX::accept(AcyclicVisitor& v) {
		Visitor<LocalCorrSurfaceABFFX>* v1 =
			dynamic_cast<Visitor<LocalCorrSurfaceABFFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABF::accept(v);
	}


	
}

