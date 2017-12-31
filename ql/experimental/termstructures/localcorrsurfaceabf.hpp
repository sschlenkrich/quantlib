/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017 Cord Harms

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

#ifndef quantlib_localcorrsurfaceabf_hpp
#define quantlib_localcorrsurfaceabf_hpp

#include <ql/experimental/termstructures/localcorrtermstructure.hpp>
#include <ql/math/interpolations/interpolation2d.hpp>

namespace QuantLib {

    //! Local Correlation surface derived 
    /*! For details about this implementation refer to
        
        \bug this class is untested, probably unreliable.
    */
    class LocalCorrSurfaceABF : public LocalCorrTermStructure {
      public:
        LocalCorrSurfaceABF(const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes,
						    const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&		  	    processToCal);
        
        
        //@}
        //! \name Visitability
        //@{
        virtual void accept(AcyclicVisitor&);
        
		//@}
		//@}
		//! \name Modifiers
		//@{
		template <class Interpolator>
		void setInterpolation(const Interpolator& i = Interpolator()) {
			interpolatorF_ =
				i.interpolate(times_.begin(), times_.end(),
					strikes_.begin(), strikes_.end(),
					surfaceF_);
			notifyObservers();
		}
		//@}
		//! \name VolatilityTermStructure interface
		//@{
		Real minStrike() const {
			return strikes_.front();
		}
		Real maxStrike() const {
			return strikes_.back();
		}
		Date maxDate() const {
			return maxDate_;
		}
      protected:
		  void localCorrImpl(RealStochasticProcess::MatA& corrMatrix, Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false) const;
		  virtual QuantLib::Real localA(Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false) const = 0;
		  virtual QuantLib::Real localB(Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false) const = 0;
		  virtual QuantLib::Real localF(Time t, const RealStochasticProcess::VecA& X0,
			  bool extrapolate = false) const = 0;

		  virtual void initializeF() =0;

		  RealStochasticProcess::MatA corr0_;
		  RealStochasticProcess::MatA corr1_;

		  Interpolation2D interpolatorF_;

		  Date maxDate_;
		  std::vector<Real> strikes_;
		  std::vector<Time> times_;
		  Matrix surfaceF_;

	  private:
        
    };

}

#endif
