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

/*! \file ParticleMethodUtils.hpp
    \brief particle method for local correlation derived ....
*/

#ifndef quantlib_particleMethodUtils_hpp
#define quantlib_particleMethodUtils_hpp

#include <ql\experimental\termstructures\Helper\KernelInterface.hpp>
#include <ql\experimental\termstructures\Helper\CalibratorLocalCorrInt.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql\handle.hpp>
#include <string.h>

namespace QuantLib {

    //!  
    /*! 
	*/
    class ParticleMethodUtils : public CalibratorLocalCorrInt {
      public:
		  //ParticleMethodUtils() : CalibratorLocalCorrInt(){};
		  ParticleMethodUtils(const std::string& kernel, unsigned int numberOfPaths, Time maxTime,
			  Time deltaT, Time tMin, Real kappa, Real sigmaAVR, Real exponentN, Real gridMinQuantile,
			  Real gridMaxQuantile, unsigned int ns1, unsigned int ns2);

		  virtual void calibrateFX(std::vector<std::vector<Real>>& strikes, std::vector<Time>& times, std::vector<std::vector<Real>>& surfaceF,
			  const std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>>& processes,
			  const boost::shared_ptr<GeneralizedBlackScholesProcess>& processToCal);
		  
		  /*boost::shared_ptr<KernelInterface>& getKernel() { return kernel_; };
		  unsigned int getNumberOfPaths() { return numberOfPaths_; };
		  Time getMaxTime() { return maxTime_;};
		  Time getDeltaT() { return deltaT_;};
		  Time getTMin() { return tMin_;};
		  Real getKappa() { return kappa_; };
		  Real getSigmaAVR() { return sigmaAVR_; };
		  Real getExponentN() { return exponentN_; };
		  Real getGridMinQuantile() { return gridMinQuantile_; };
		  Real getGridMaxQuantile() { return gridMaxQuantile_; };*/
      protected:
		  
	  private:
		  Real bandwidth(Time t, Real s0) const;
		  Real kernel(Real bandwidth, Real x) const;
		  
		  boost::shared_ptr<KernelInterface> kernel_;
		  unsigned int numberOfPaths_;
		  Time maxTime_;
		  Time deltaT_;
		  Time tMin_;
		  Real kappa_;
		  Real sigmaAVR_;
		  Real exponentN_;
		  Real gridMinQuantile_;
		  Real gridMaxQuantile_;
		  unsigned int ns1_;
		  unsigned int ns2_;


    };

}

#endif
