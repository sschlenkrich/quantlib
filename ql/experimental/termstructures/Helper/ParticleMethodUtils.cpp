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

/*! \file localcorrsurface.hpp
    \brief Local Correlation surface derived ....
*/

#include <ql\experimental\termstructures\Helper\ParticleMethodUtils.hpp>

namespace QuantLib {

	ParticleMethodUtils::ParticleMethodUtils(const std::string& kernel, unsigned int numberOfPaths, Time maxTime,
		Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
		Real gridMaxQuantile) :
		numberOfPaths_(numberOfPaths), maxTime_(maxTime), deltaT_(deltaT), tMin_(tMin), kappa_(kappa),
		sigmaAVR_(sigmaAVR), exponentN_(exponentN), gridMinQuantile_(gridMinQuantile), gridMaxQuantile_(gridMaxQuantile) {
		
		if (kernel == "QuarticKernel") {
			kernel_ = boost::shared_ptr<KernelInterface>(new QuarticKernel());
		}
		else {
			QL_REQUIRE(false, "Kernel not supported. Supported is: QuarticKernel");
		}
	}

	/*void ParticleMethodUtils::calibrateFX(std::vector<Real>& strikes, std::vector<Time>& times, Matrix& surfaceF,
		const std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<GeneralizedBlackScholesProcess>& processToCal) {
	
	}*/	  
}
