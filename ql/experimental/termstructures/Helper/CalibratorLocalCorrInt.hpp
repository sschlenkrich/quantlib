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

/*! \file CalibratorLocalCorrInt.hpp
*/

#ifndef quantlib_calibratorlocalcorrint_hpp
#define quantlib_calibratorlocalcorrint_hpp

#include <vector>
#include <ql\math\matrix.hpp>
#include <ql\types.hpp>
#include <boost\smart_ptr\shared_ptr.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    //! CalibratorLocalCorrInt derived 
    /*! For details about this implementation refer to
        
    */
    class CalibratorLocalCorrInt {
      public:
		  virtual void calibrateFX(std::vector<std::vector<Real>>& strikes, std::vector<Time>& times, std::vector<std::vector<Real>>& surfaceF,
			  const std::vector<boost::shared_ptr<GeneralizedBlackScholesProcess>>& processes,
			  const boost::shared_ptr<GeneralizedBlackScholesProcess>& processToCal) = 0;

      protected:
		  
	  private:
			
    };

}

#endif
