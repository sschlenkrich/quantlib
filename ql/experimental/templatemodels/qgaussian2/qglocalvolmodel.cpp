/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018, Sebastian Schlenkrich

*/

#include <ql/pricingengines/blackformula.hpp>
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
		const BigNatural                                       seed,
		const size_t                                           debugLevel)
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
		volTS_(volTS), swapIndex_(swapIndex), stdDevGrid_(stdDevGrid), nPaths_(nPaths), seed_(seed), debugLevel_(debugLevel){
        // we can't calibrate here because *this needs to be assigned to shared pointer first
	}

	// we cache the float leg and fixed leg valuation to speed-up MC simulation
	QGLocalvolModel::SwaptionFactory::SwaptionFactory(const Time obsTime, const SwapCashFlows& scf)
		: floatLeg_(boost::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCPayoff::Cache(
			boost::shared_ptr<QGLocalvolModel::MCPayoff>(new MCAnnuity(obsTime, scf.floatTimes(), scf.floatWeights()))))),
		annuityLeg_(boost::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCPayoff::Cache(
			boost::shared_ptr<QGLocalvolModel::MCPayoff>(new MCAnnuity(obsTime, scf.fixedTimes(), scf.annuityWeights())))))
	{}

	boost::shared_ptr<QGLocalvolModel::MCPayoff> QGLocalvolModel::SwaptionFactory::swaption(const Real strike, const Real callOrPut) {
		boost::shared_ptr<QGLocalvolModel::MCPayoff> fixedRate(new QGLocalvolModel::MCPayoff::FixedAmount(strike));
		boost::shared_ptr<QGLocalvolModel::MCPayoff> fixedLeg(new QGLocalvolModel::MCPayoff::Mult(fixedRate, annuityLeg_));
		boost::shared_ptr<QGLocalvolModel::MCPayoff> swap;
		if (callOrPut == 1.0) {
			swap = boost::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCPayoff::Axpy(-1.0, fixedLeg, floatLeg_));
		}
		else {
			swap = boost::shared_ptr<QGLocalvolModel::MCPayoff>(new QGLocalvolModel::MCPayoff::Axpy(-1.0, floatLeg_, fixedLeg));
		}
		boost::shared_ptr<QGLocalvolModel::MCPayoff> zero(new QGLocalvolModel::MCPayoff::FixedAmount(0.0));
		boost::shared_ptr<QGLocalvolModel::MCPayoff> swpt(new QGLocalvolModel::MCPayoff::Max(swap, zero));
		boost::shared_ptr<QGLocalvolModel::MCPayoff> pay(new QGLocalvolModel::MCPayoff::Pay(swpt, floatLeg_->observationTime()));
		return pay;
	}

	void QGLocalvolModel::simulateAndCalibrate() {
		debugLog_.clear();
		// reset local volatility
		sigmaS_.clear();
		strikeGrid_.clear();
		locvolGrid_.clear();
		// initialise MC simulation
		simulation_ = boost::shared_ptr<MCSimulation>(new MCSimulation(shared_from_this(), this->times(), this->times(), nPaths_, seed_, false, true, true));
		simulation_->prepareSimulation();
		simulation_->simulate(0);
		// prepare for first simulation step using approximate local vol
		boost::shared_ptr<QuasiGaussianModel> thisAsQGModel = boost::static_pointer_cast<QuasiGaussianModel>(shared_from_this());
		Date today = termStructure()->referenceDate(); // check if this is the correct date...
		// we need to set up an initial swap rate model for gradient calculation
		Date fixingDate = swapIndex_->fixingCalendar().advance(today, 1, Days, Following);
		SwapCashFlows scf(swapIndex_->underlyingSwap(fixingDate), termStructure_, true); 
		std::vector<Real> swapRateModelTimes{ 0.0, Actual365Fixed().yearFraction(today,fixingDate) };
		sigmaMode_ = Parent; // switch-off sigma_x calculation; we only need yield curve information
		swapRateModel_ = boost::shared_ptr<QGSwaprateModel>(new QGSwaprateModel(thisAsQGModel,
			scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights(), swapRateModelTimes, false));
		sigmaMode_ = Calibration; // switch-on sigma_x calculation for simulation and evolve calls
		// run simulation and calibration
		for (size_t idx = 0; idx < times().size(); ++idx) {
			// simulate next step T[idx-1] to T[idx], note X[] starts at 0
			simulation_->simulate(idx+1);

			// specify swap rate at T[idx]
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
			Real stdDev = smileSection->optionPrice(swapRate, Option::Call) / M_1_SQRTPI / M_SQRT_2;
			std::vector<Real> initialStrikeGrid(stdDevGrid_);
			for (size_t k = 0; k < initialStrikeGrid.size(); ++k) initialStrikeGrid[k] = initialStrikeGrid[k] * stdDev + swapRate;

			Time obsTime = (idx>0) ? (times()[idx-1]) : (0.0);
			SwaptionFactory calibFactory(obsTime, scf);
			SwaptionFactory testFactory(times()[idx], scf);
			// include strike-adjustment to ensure put-call parity in simulation...
			boost::shared_ptr<MCPayoff> mcAtmCall = calibFactory.swaption(swapRate, 1.0);
			boost::shared_ptr<MCPayoff> mcAtmPut  = calibFactory.swaption(swapRate, -1.0);
			Real atmCall = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcAtmCall), simulation_);
			Real atmPut  = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcAtmPut),  simulation_);
			Real nu = (atmPut - atmCall) / annuity;
			// logging
			if (debugLevel_>0) debugLog_.push_back("T = " + std::to_string(times()[idx]) + ", fixingDate = " + std::to_string(fixingDate.serialNumber()) + ", swapRate = " + std::to_string(swapRate) + ", annuity = " + std::to_string(annuity) + ", stdDev = " + std::to_string(stdDev) + ", nu = " + std::to_string(nu));

			// these are the actual grids used for local vol calculation
			std::vector<Real> strikeGrid;
			std::vector<Real> localVol;
			for (size_t k = 0; k < initialStrikeGrid.size(); ++k) {

				Real h = 1.0e-5; // 0.1bp shift size
				Real call = smileSection->optionPrice(initialStrikeGrid[k], Option::Call);
				Real callph = smileSection->optionPrice(initialStrikeGrid[k] + h, Option::Call);
				Real callmh = smileSection->optionPrice(initialStrikeGrid[k] - h, Option::Call);

				// we only calculate local vol if within [1e-6, 1 - 1e-6] quantile 
				Real quantile = 1.0 + (callph - callmh) / 2.0 / h;
				if (debugLevel_>1) debugLog_.push_back("k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", call = " + std::to_string(call) + ", callph = " + std::to_string(callph) + ", callmh = " + std::to_string(callmh) + ", quantile = " + std::to_string(quantile));
				if ((quantile < 1.0e-6) || (quantile > 1.0 - 1.e-6)) {
					if (debugLevel_>1) debugLog_.push_back("Warning: Skip local vol.");
					continue;  // skip calculation
				}

				boost::shared_ptr<MCPayoff> mcCall = calibFactory.swaption(initialStrikeGrid[k] - nu, 1.0);
				Real mcCallNpv = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcCall), simulation_);
				mcCallNpv /= annuity;
				if (debugLevel_ > 1) {
					boost::shared_ptr<MCPayoff> mcPut = calibFactory.swaption(initialStrikeGrid[k] - nu, -1.0);
					Real mcPutNpv = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcPut), simulation_);
					mcPutNpv /= annuity;
					debugLog_.push_back("k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", mcCallNpv = " + std::to_string(mcCallNpv) + ", mcPutNpv = " + std::to_string(mcPutNpv) + ", C/P-parity = " + std::to_string(mcCallNpv - mcPutNpv - (swapRate - initialStrikeGrid[k])));
				}

				Real d2CalldK2 = (callph + callmh - 2.0*call) / h / h;
				if (d2CalldK2 < 1.0e-8) {
					if (debugLevel_>1) debugLog_.push_back("Warning: d2CalldK2 = " + std::to_string(d2CalldK2) + "Skip local vol.");
					continue;  // skip calculation
				}
				Real dCalldT = (call - mcCallNpv) / (times()[idx] - obsTime);
				if (dCalldT < 1.0e-8) {
					if (debugLevel_>1) debugLog_.push_back("Warning: dCalldT = " + std::to_string(dCalldT) + "Skip local vol.");
					continue;  // skip calculation
				}
				Real lvol    = sqrt(2 * dCalldT / d2CalldK2);
				if (debugLevel_>1) debugLog_.push_back("k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", d2CalldK2 = " + std::to_string(d2CalldK2) + ", dCalldT = " + std::to_string(dCalldT) + ", lvol = " + std::to_string(lvol));
				strikeGrid.push_back(initialStrikeGrid[k]);
				localVol.push_back(lvol);

				// test calibration
				if (debugLevel_ > 1) {
					boost::shared_ptr<MCPayoff> mcCall = testFactory.swaption(initialStrikeGrid[k], 1.0);
					boost::shared_ptr<MCPayoff> mcPut = testFactory.swaption(initialStrikeGrid[k], -1.0);
					Real testCall = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcCall), simulation_);
					Real testPut = MCPayoff::Pricer::NPV(std::vector< boost::shared_ptr<MCPayoff> >(1, mcPut), simulation_);
					testCall /= annuity;
					testPut /= annuity;
					debugLog_.push_back("k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", testCall = " + std::to_string(testCall) + ", testPut = " + std::to_string(testPut) + ", call = " + std::to_string(call));
					try {
						Real callVol = bachelierBlackFormulaImpliedVol(Option::Call, initialStrikeGrid[k], swapRate, times()[idx], testCall);
						Real putVol = bachelierBlackFormulaImpliedVol(Option::Put, initialStrikeGrid[k], swapRate, times()[idx], testPut);
						Real vanVol = bachelierBlackFormulaImpliedVol(Option::Call, initialStrikeGrid[k], swapRate, times()[idx], call);
						debugLog_.push_back("k = " + std::to_string(k) + ", strike = " + std::to_string(initialStrikeGrid[k]) + ", callVol = " + std::to_string(callVol) + ", putVol = " + std::to_string(putVol) + ", vanVol = " + std::to_string(vanVol));
					}
					catch (std::exception e) {
						std::string what = e.what();
						debugLog_.push_back("Error: " + what);
					}
				}
				// done!
			}

			// set up interpolation
			// sigmaS_.push_back(LinearInterpolation(strikeGrid.begin(), strikeGrid.end(), localVol.begin()));
			strikeGrid_.push_back(strikeGrid);
			locvolGrid_.push_back(localVol);
			sigmaS_.push_back(Linear().interpolate(strikeGrid_.back().begin(), strikeGrid_.back().end(), locvolGrid_.back().begin()));

			// set up swap rate model here such that we don't need to care about it during sigma_x calculation
			std::vector<Real> swapRateModelTimes{ 0.0, times()[idx] };  // gradient observation date is at next grid point
			sigmaMode_ = Parent; // switch-off sigma_x calculation; we only need yield curve information
			swapRateModel_ = boost::shared_ptr<QGSwaprateModel>(new QGSwaprateModel(thisAsQGModel,
				scf.floatTimes(), scf.floatWeights(), scf.fixedTimes(), scf.annuityWeights(), swapRateModelTimes, false));
			sigmaMode_ = Calibration; // switch-on sigma_x calculation for simulation and evolve calls
		}
	}

	inline std::vector< std::vector<Real> >
	QGLocalvolModel::sigma_xT(const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >&  y) {
		if (sigmaMode_ == Parent) return QuasiGaussianModel::sigma_xT(t, x, y);  // this is more like a fall back if volatility is irrelevant
		if (sigmaMode_ == Calibration) {
			// the current implementation is intended only for calibration phase
			// we rely on the swap rate model beeing initialised and updated properly before this function is called
			// this will most likely not work if model is used with external MC simulation
			QL_REQUIRE(swapRateModel_ != 0, "swapRateModel_!=0 required");
			Real observationTime = swapRateModel_->modelTimes().back();  // this should be times()[idx] from calibrateAndSimulate()
			Real swapRate = swapRateModel_->swapRate(observationTime, x, y);
			Real swapGradient = (swapRateModel_->swapGradient(observationTime, x, y))[0];
			Real lvol;
			if (sigmaS_.size() == 0) {
				lvol = volTS_->volatility(observationTime, swapIndex_->tenor(), swapRate); // initial approximation
			}
			else {  // ensure flat extrapolation..
				bool warning = false;
				if (swapRate < sigmaS_.back().xMin()) {
					if (debugLevel_ > 2) debugLog_.push_back("Warning: xMin = " + std::to_string(sigmaS_.back().xMin()) + ", swapRate = " + std::to_string(swapRate));
					swapRate = sigmaS_.back().xMin();
					warning = true;
				}
				if (swapRate > sigmaS_.back().xMax()) {
					if (debugLevel_ > 2) debugLog_.push_back("Warning: xMax = " + std::to_string(sigmaS_.back().xMax()) + ", swapRate = " + std::to_string(swapRate));
					swapRate = sigmaS_.back().xMax();
					warning = true;
				}
				lvol = sigmaS_.back()(swapRate);  // this relies on the usage within calibrateAndSimulate()													   
				if ((warning)&&(debugLevel_ > 2)) {
					debugLog_.push_back("Warning: lvol = " + std::to_string(lvol) + ", swapGradient = " + std::to_string(swapGradient));
				}
			}
			if (debugLevel_ > 3) debugLog_.push_back("t = " + std::to_string(t) +
			                                         ", x = " + std::to_string(x[0]) +
				                                     ", y = " + std::to_string(y[0][0]) +
				                                     ", swapRate = " + std::to_string(swapRate) +
				                                     ", swapGradient = " + std::to_string(swapGradient) +
				                                     ", lvol = " + std::to_string(lvol) );
			return  std::vector< std::vector<Real> >(1, std::vector<Real>(1, lvol / swapGradient));
		}
		QL_REQUIRE(false, "post-calibration mode for sigma_x not implemented yet");
	}

}

