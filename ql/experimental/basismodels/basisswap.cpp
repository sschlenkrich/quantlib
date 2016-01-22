/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswap.hpp
    \brief (Cross currency) Interest rate swap consisting of two swaps
*/


#include <ql/experimental/basismodels/basisswap.hpp>



namespace QuantLib {
    
	Real BasisSwap::fairRate() const {
		Real basisPoint = 1.0e-4;
		Real num = ( (calcParSpread_) ? 0.0 : legNPV(parLegIndex_) ) - NPV();
		Real den = legBPS(parLegIndex_) / basisPoint;
		return num / den;
	}

}


