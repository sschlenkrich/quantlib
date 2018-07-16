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
		typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::Annuity MCAnnuity;
		typedef RatesPayoffT<QuantLib::Time, QuantLib::Real, QuantLib::Real>::GeneralSwaption MCSwaption;

	private:
		// calibration target volatilities
		boost::shared_ptr<SwaptionVolatilityStructure> volTS_;
        
		// simulation derived during calibration; we need to store nPaths_ and seed_ separately
		// because we can't initialise simulation within constructor
		size_t     nPaths_;
		BigNatural seed_;
		boost::shared_ptr<MCSimulation> simulation_;
        
		// we calibrate to a strip of swaption volatilities; maybe also co-terminals can be relevant
		boost::shared_ptr<SwapIndex> swapIndex_;

        // we describe local volatility sigmaS as a set of 1-D interpolations per time step 1 to N (excluding 0)
        std::vector< Interpolation > sigmaS_;
		std::vector < std::vector<Real> > strikeGrid_;  // we need to separately store data for interpolation
		std::vector < std::vector<Real> > locvolGrid_;


		// we have three modes for sigma_x calculation (during and after calibration phase)
		enum { 
			Parent,         // delegate sigma_x calculation to parent class; this typically doesn't do something meaningfull
			Calibration,    // calculate sigma_x as specified for calibration procedure
			Pricing         // use calibrated local vol after it is fully calibrated
		} sigmaMode_;

		// we specify the local vol grid in terms of standard deviations
		std::vector<Real> stdDevGrid_;

        // we need to cache the swap rate model and observation time for the current time step
		boost::shared_ptr<QGSwaprateModel> swapRateModel_;

		// we need a factory for MC swaptions for efficient calculation
		class SwaptionFactory {
			boost::shared_ptr<MCPayoff> floatLeg_;
			boost::shared_ptr<MCPayoff> annuityLeg_;
		public:
			SwaptionFactory(const Time obsTime, const SwapCashFlows& scf);
			boost::shared_ptr<MCPayoff> swaption(const Real strike, const Real CallOrPut);
		};

		// debugging, warning and errors
		std::vector<std::string> debugLog_;
		size_t debugLevel_; // 0 ... no debugging
		                    // 1 ... debugging for time steps
		                    // 2 ... debugging for strikes
		                    // 3 ... warnings for simulation
		                    // 4 ... debugging for each path (not recommended)

	public:

		QGLocalvolModel( const Handle<YieldTermStructure>&                      termStructure,
			             const boost::shared_ptr<SwaptionVolatilityStructure>&  volTS,
			             const Real                                             chi,
			             const boost::shared_ptr<SwapIndex>&                    swapIndex,
			             const std::vector<Real>&                               times,
			             const std::vector<Real>&                               stdDevGrid,
			             const size_t                                           nPaths,
			             const BigNatural                                       seed = 1234,
			             const size_t                                           debugLevel = 1);

		// do the actual calculation
		void simulateAndCalibrate();

		inline virtual std::vector< std::vector<Real> >
		sigma_xT(const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >&  y);

		// inspectors
		inline const boost::shared_ptr<MCSimulation> simulation() { return simulation_; }
		inline const Real sigmaS(const size_t idx, const Real s) { return sigmaS_[idx](s); }

		std::vector<std::string> debugLog() { return debugLog_; }

		// test the calibration of the model
		std::vector< std::vector<Real> > calibrationTest(const std::vector<Date>&  exerciseDates,
			                                             const std::vector<Real>&  stdDevStrikes );
    };

}

#endif  /* ifndef quantlib_qglocalvolmodel_hpp */
