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

#include <ql\experimental\termstructures\localCorrFX\CTSlocalInCrossCorrelationFX.hpp>
#include <ql/math/interpolations/bilinearinterpolation.hpp>

#include <ql\experimental\termstructures\Helper\ParticleMethodUtils.hpp>

namespace QuantLib {

	CTSlocalInCrossCorrelationFX::CTSlocalInCrossCorrelationFX(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal,
		boost::shared_ptr<CalibratorLocalCorrInt>&										calibratorLocalCorr)
    : LocalCorrSurfaceABFFX(processes,processToCal,calibratorLocalCorr){
		initializeF();
		setInterpolation<Bilinear>();

		//ParticleMethodUtils test = ParticleMethodUtils(std::string(""),2,Time(2), Time(2), Time(2),Real(2), Real(2), Real(2), Real(2), Real(2));
    }

	QuantLib::Real CTSlocalInCrossCorrelationFX::localA(Time t, const RealStochasticProcess::VecA& X0,
		bool extrapolate) const {
		return 0;
	}

	QuantLib::Real CTSlocalInCrossCorrelationFX::localB(Time t, const RealStochasticProcess::VecA& X0,
		bool extrapolate) const {
		return 1;
	}
	void CTSlocalInCrossCorrelationFX::accept(AcyclicVisitor& v) {
		Visitor<CTSlocalInCrossCorrelationFX>* v1 =
			dynamic_cast<Visitor<CTSlocalInCrossCorrelationFX>*>(&v);
		if (v1 != 0)
			v1->visit(*this);
		else
			LocalCorrSurfaceABFFX::accept(v);
	}


	void CTSlocalInCrossCorrelationFX::initializeF() {

		//calibratorLocalCorr_->calibrateFX(strikes_, times_, surfaceF_, processes_, processToCal_);

		//strikes_.resize(2);
		//times_.resize(2);
		//surfaceF_ = Matrix(2,2);

		//for (size_t i = 0; i < surfaceF_.size1(); i++)
		//{	
		//	for (size_t j = 0; j < surfaceF_.size2(); j++)
		//	{
		//		surfaceF_[i][j] = 0.9;
		//	}
		//}

		//times_[0] = 0;
		//times_[1] = 1;
		//strikes_[0] = 0.7;
		//strikes_[1] = 1;
	}
}

