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

#ifndef quantlib_kernelInterface_hpp
#define quantlib_kernelInterface_hpp


namespace QuantLib {

    /**/
    class KernelInterface {
      public:
		  KernelInterface() {};
		  virtual Real operator()(Real x) const = 0;
      protected:
		  
	  private:
        
    };

	class QuarticKernel : public KernelInterface {
	public:
		QuarticKernel() : KernelInterface() {};
		virtual Real operator()(Real x) const { return abs(x) < 1 ? (1 - x*x) : 0; };
	};

}

#endif
