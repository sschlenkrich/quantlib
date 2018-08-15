/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018, Sebastian Schlenkrich

*/



#ifndef quantlib_qglsvmodel_hpp
#define quantlib_qglsvmodel_hpp

#include <ql/experimental/templatemodels/qgaussian2/qglocalvolmodel.hpp>


namespace QuantLib {


	// alternative backward-looking calibration methodology for QGLocalvolModel
	class QGLSVModel : public QGLocalvolModel {
    protected:
        size_t nStrikes_;
	public:
		QGLSVModel(
			const Handle<YieldTermStructure>&                      termStructure,
			const boost::shared_ptr<SwaptionVolatilityStructure>&  volTS,
			const Real                                             chi,
			const Real                                             theta,
			const Real                                             eta,
			const boost::shared_ptr<SwapIndex>&                    swapIndex,
			const std::vector<Real>&                               times,
			const size_t                                           nStrikes,
			const bool                                             calcStochVolAdjustment,
			const size_t                                           nPaths,
			const BigNatural                                       seed = 1234,
			const size_t                                           debugLevel = 1)
			: QGLocalvolModel(termStructure, volTS, chi, theta, eta, swapIndex, times, std::vector<Real>(), calcStochVolAdjustment, 0.0, nPaths, seed, debugLevel), nStrikes_(nStrikes) {}
		// do the actual calculation
        
		virtual void simulateAndCalibrate();
	};




}

#endif  /* ifndef quantlib_qglsvmodel_hpp */
