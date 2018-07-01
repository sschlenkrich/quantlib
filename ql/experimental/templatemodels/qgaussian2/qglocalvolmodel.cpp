/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018, Sebastian Schlenkrich

*/


#include <ql/math/interpolations/linearinterpolation.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>

#include <ql/experimental/templatemodels/qgaussian2/qglocalvolmodel.hpp>

namespace QuantLib {

    // model implementaion
	QGLocalvolModel::QGLocalvolModel(
		const Handle<YieldTermStructure>&                      termStructure,
		const boost::shared_ptr<SwaptionVolatilityStructure>&  volTS,
		const Real                                             chi,
		const boost::shared_ptr<SwapIndex>&                    swapIndex,
		const std::vector< Real >&                             times,
		const std::vector<Real>&                               stdDevGrid,
		const size_t                                           nPaths,
		const BigNatural                                       seed)
		: QuasiGaussianModel(termStructure, 1, times,
			std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // sigma
			std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // slope
			std::vector< std::vector<Real> >(1, std::vector<Real>(times.size(), 0.0)),  // curve
			std::vector<Real>(times.size(), 0.0),                                       // eta
			std::vector<Real>(1, 1.0),                                                  // delta
			std::vector<Real>(1, chi),                                                  // mean reversion
			std::vector< std::vector<Real> >(1, std::vector<Real>(1, 1.0)),             // Gamma
			0.1                                                                         // theta
		),
		volTS_(volTS), swapIndex_(swapIndex), stdDevGrid_(stdDevGrid), sigmaSIsCalibrated_(false) {
		simulation_ = boost::shared_ptr<MCSimulation>(new MCSimulation(shared_from_this(), this->times(), this->times(), nPaths, seed, true, true, false));
		calibrateAndSimulate();
	}

	inline void QGLocalvolModel::calibrateAndSimulate() {
		simulation_->prepareSimulation();
		simulation_->simulate(0);
		for (size_t idx = 1; idx < times().size(); ++idx) {
			// specify swap rate
			Date today = termStructure()->referenceDate(); // check if this is the correct date...
			Date fixingDate = today + ((BigInteger)ClosestRounding(0)(times()[idx]*365.0)); // assuming act/365 day counting
			if (!swapIndex_->isValidFixingDate(fixingDate))	fixingDate = swapIndex_->fixingCalendar().adjust(fixingDate, Following);
			SwapCashFlows scf(swapIndex_->underlyingSwap(fixingDate), termStructure_, true);        // assume continuous tenor spreads
			Real annuity = 0.0;
			for (size_t k = 0; k < scf.fixedTimes().size(); ++k) annuity += scf.annuityWeights()[k] * termStructure()->discount(scf.fixedTimes()[k]);
			Real floatLeg = 0.0;
			for (size_t k = 0; k < scf.floatTimes().size(); ++k) floatLeg += scf.floatWeights()[k] * termStructure()->discount(scf.floatTimes()[k]);
			Real swapRate = floatLeg / annuity;

			// set up smile section and strike grid
			boost::shared_ptr<SmileSection> smileSection = volTS_->smileSection(times()[idx], swapIndex_->tenor(), true);
			Real atmCall = smileSection->optionPrice(swapRate, Option::Call);
			Real stdDev = atmCall / M_1_SQRTPI / M_SQRT_2;
			std::vector<Real> strikeGrid(stdDevGrid_);
			for (size_t k = 0; k < strikeGrid.size(); ++k) strikeGrid[k] = strikeGrid[k] / stdDev + swapRate;

			// calculate call prices; this can be streamlined, no need to store all vectors
			std::vector< boost::shared_ptr<MCPayoff> > mcSwaptions(strikeGrid.size());
			std::vector<Real> npvMC(strikeGrid.size());
			std::vector<Real> localVol(strikeGrid.size());
			for (size_t k = 0; k < strikeGrid.size(); ++k) {
				mcSwaptions[k] = boost::shared_ptr<MCPayoff>(new MCSwaption(times()[idx-1], scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights(), strikeGrid[k], 1.0));
		        npvMC[k] = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcSwaptions[k]), simulation_);
				npvMC[k] /= annuity;
				Real h = 1.0e-5; // 0.1bp shift size
				Real call   = smileSection->optionPrice(strikeGrid[k], Option::Call);
				Real callph = smileSection->optionPrice(strikeGrid[k]+h, Option::Call);
				Real callmh = smileSection->optionPrice(strikeGrid[k]-h, Option::Call);
				Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
				Real dCalldT = (call - npvMC[k]) / (times()[idx] - times()[idx-1]);
				localVol[k] = sqrt(2 * dCalldT / d2CalldK2);
				// ensure local vol is positive...
			}

			// set up interpolation
			sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));

			// set up swap rate model here such that we don't need to care about it during sigma_x calculation
			boost::shared_ptr<QuasiGaussianModel> thisAsQGModel = boost::static_pointer_cast<QuasiGaussianModel>(shared_from_this());
			std::vector<Real> swapRateModelTimes{ 0.0, times()[idx] };  // gradient observation date is at next grid point
			swapRateModel_ = boost::shared_ptr<QGSwaprateModel>(new QGSwaprateModel(thisAsQGModel,
				scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights(), swapRateModelTimes, false));

			// simulate next step T[idx-1] to T[idx]
			simulation_->simulate(idx);
		}
	}

	inline std::vector< std::vector<Real> >
	QGLocalvolModel::sigma_xT(const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >&  y) {
		// the current implementation is intended only for calibration phase
		QL_REQUIRE(!sigmaSIsCalibrated_, "post-calibration mode for sigma_x not implemented yet");
		// we rely on the swap rate model beeing initialised and updated properly before this function is called
		// this will most likely not work if model is used with external MC simulation
		QL_REQUIRE(swapRateModel_ != 0, "swapRateModel_!=0 required");
		Real observationTime = swapRateModel_->modelTimes().back();  // this should be times()[idx] from calibrateAndSimulate()
		Real swapRate = swapRateModel_->swapRate(observationTime, x, y);
		Real swapGradient = (swapRateModel_->swapGradient(observationTime, x, y))[0];
		Real lvol = sigmaS_.back()(swapRate);  // this relies on the usage within calibrateAndSimulate()
		// ensure proper extrapolation....
		return  std::vector< std::vector<Real> >(1, std::vector<Real>(1, lvol / swapGradient));
	}

}

