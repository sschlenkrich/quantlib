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

namespace QuantLib {

	CTSlocalInCrossCorrelationFX::CTSlocalInCrossCorrelationFX(
		const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&			    processToCal)
    : LocalCorrSurfaceABFFX(processes,processToCal){
		
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
}

