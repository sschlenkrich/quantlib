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
#include <math.h>

namespace QuantLib {

	ParticleMethodUtils::ParticleMethodUtils(const std::string& kernel, unsigned int numberOfPaths, Time maxTime,
		Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
		Real gridMaxQuantile, unsigned int ns1, unsigned int ns2) :
		numberOfPaths_(numberOfPaths), maxTime_(maxTime), deltaT_(deltaT), tMin_(tMin), kappa_(kappa),
		sigmaAVR_(sigmaAVR), exponentN_(exponentN), gridMinQuantile_(gridMinQuantile), 
		gridMaxQuantile_(gridMaxQuantile), ns1_(ns1), ns2_(ns2) {
	
		if (kernel == "QuarticKernel") {
			kernel_ = boost::shared_ptr<KernelInterface>(new QuarticKernel());
		}
		else {
			QL_REQUIRE(false, "Kernel not supported. Supported is: QuarticKernel");
		}
	}

	void ParticleMethodUtils::calibrateFX(std::vector<std::vector<Real>>& strikes, std::vector<Time>& times, std::vector<std::vector<Real>>& surfaceF,
		const std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>>& processes,
		const boost::shared_ptr<GeneralizedBlackScholesProcess>& processToCal) {
		
		//time grid from t=0 to maxTime:
		times.resize(1);
		times[0] = 0;
		size_t i = 0;
		while (times[i] * deltaT_ < maxTime_) {
			times.push_back(times[i]+deltaT_);
			i++;
		}

		//start to create strike grid

		surfaceF.resize(times.size());
		strikes.resize(times.size());

		for (size_t i = 0; i < surfaceF.size(); i++)
		{
			surfaceF[i].resize(2);
			strikes[i].resize(2);
			for (size_t j = 0; j < surfaceF[i].size(); j++)
			{
				surfaceF[i][j] = 0.9;
				strikes[i][0] = 0.7;
				strikes[i][1] = 1.2;
			}
		}
	}	  

	Real ParticleMethodUtils::bandwidth(Time t, Real s0) const {
		return kappa_*sigmaAVR_* s0 * sqrt(t>tMin_ ? t : tMin_) * pow(numberOfPaths_,exponentN_);
	}

	Real ParticleMethodUtils::kernel(Real bandwidth, Real x) const {
		QL_REQUIRE(bandwidth != 0, "Error in ParticleMethodUtils: bandwidth is not allowed to be zero.");
		return kernel_->value(x / bandwidth) / bandwidth;
	}
}
