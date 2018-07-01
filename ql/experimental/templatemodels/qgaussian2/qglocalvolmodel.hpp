/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018, Sebastian Schlenkrich

*/



#ifndef quantlib_qglocalvolmodel_hpp
#define quantlib_qglocalvolmodel_hpp

#include <ql/math/interpolation.hpp>

#include <ql/indexes/swapindex.hpp>
#include <ql/instruments/swaption.hpp>
#include <ql/termstructures/volatility/swaption/swaptionvolstructure.hpp>

//#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>

#include <ql/experimental/basismodels/swaptioncfs.hpp>

#include <ql/experimental/templatemodels/qgaussian2/quasigaussianmodel2T.hpp>
#include <ql/experimental/templatemodels/qgaussian2/qgswapratemodelT.hpp>
#include <ql/experimental/templatemodels/montecarlo/mcsimulationT.hpp>
#include <ql/experimental/templatemodels/montecarlo/ratespayoffT.hpp>


namespace QuantLib {

	// calibrate Quasi Gaussian model to implied normal volatilities
    class QGLocalvolModel : public QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> { 
	public:
	    typedef QuasiGaussianModel2T<QuantLib::Time, QuantLib::Real, QuantLib::Real> QuasiGaussianModel;
	    typedef QGSwaprateModelT<QuantLib::Time, QuantLib::Real, QuantLib::Real> QGSwaprateModel;
	    typedef MCSimulationT<QuantLib::Time,QuantLib::Real,QuantLib::Real> MCSimulation;    
		typedef MCPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real> MCPayoff;
		typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::GeneralSwaption MCSwaption;

	private:
		// calibration target volatilities
		boost::shared_ptr<SwaptionVolatilityStructure> volTS_;
        
		// calibration target volatilities
		boost::shared_ptr<MCSimulation> simulation_;
        
		// we calibrate to a strip of swaption volatilities; maybe also co-terminals can be relevant
		boost::shared_ptr<SwapIndex> swapIndex_;

        // we describe local volatility sigmaS as a set of 1-D interpolations per time step 1 to N (excluding 0)
        std::vector< Interpolation > sigmaS_;

		// we have two modes for sigma_x calculation (during and after calibration phase)
		bool sigmaSIsCalibrated_;

		// we specify the local vol grid in terms of standard deviations
		std::vector<Real> stdDevGrid_;

        // we need to cache the swap rate model and observation time for the current time step
		boost::shared_ptr<QGSwaprateModel> swapRateModel_;

		// do the actual calculation
		inline void calibrateAndSimulate();

	public:

		QGLocalvolModel( const Handle<YieldTermStructure>&                      termStructure,
			             const boost::shared_ptr<SwaptionVolatilityStructure>&  volTS,
			             const Real                                             chi,
			             const boost::shared_ptr<SwapIndex>&                    swapIndex,
			             const std::vector<Real>&                               times,
			             const std::vector<Real>&                               stdDevGrid,
			             const size_t                                           nPaths,
			             const BigNatural                                       seed = 1234 );

		inline virtual std::vector< std::vector<Real> >
		sigma_xT(const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >&  y);

		// inspectors
		inline const boost::shared_ptr<MCSimulation> simulation() { return simulation_; }
		inline const Real sigmaS(const size_t idx, const Real s) { return sigmaS_[idx](s); }

    };

}

#endif  /* ifndef quantlib_qglocalvolmodel_hpp */
