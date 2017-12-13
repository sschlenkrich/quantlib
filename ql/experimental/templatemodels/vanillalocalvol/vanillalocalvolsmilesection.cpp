/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Sebastian Schlenkrich

*/

#include <ql/experimental/templatemodels/vanillalocalvol/vanillalocalvolsmilesection.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/pricingengines/blackformula.hpp>

#include <ql/experimental/templatemodels/auxilliaries/auxilliariesT.hpp>

namespace QuantLib {

	namespace {

		class VanillaLocalVolHelper : public CostFunction {
			boost::shared_ptr<VanillaLocalVolModel>     model_;
			std::vector<Rate>                           relativeStrikes_;
			std::vector<Real>                           smilePrices_;
			std::vector<Real>                           vanillaVega_;
			Real                                        extrapolationRelativeStrike_;
			Real                                        extrapolationSlope_;
			Real                                        minSlope_;   //  lower boundary for m in calibration
			Real                                        maxSlope_;   //  upper boundary for m in calibration
			Real                                        alpha_;      //  Tikhonov alpha
			size_t                                      zeroIdx_;    //  the index with relativeStrike[.] = 0
		public:
			VanillaLocalVolHelper(
				const boost::shared_ptr<VanillaLocalVolModel>&     model,
				const std::vector<Rate>&                           relativeStrikes,
				const std::vector<Real>&                           smilePrices,
				const std::vector<Real>&                           vanillaVega,
				const Real                                         extrapolationRelativeStrike,
				const Real                                         extrapolationSlope,
				const Real                                         minSlope,   //  lower boundary for m in calibration
				const Real                                         maxSlope,   //  upper boundary for m in calibration
				const Real                                         alpha)     //  Tikhonov alpha
				: model_(model), relativeStrikes_(relativeStrikes), smilePrices_(smilePrices), vanillaVega_(vanillaVega),
				extrapolationRelativeStrike_(extrapolationRelativeStrike), extrapolationSlope_(extrapolationSlope),
				minSlope_(minSlope), maxSlope_(maxSlope), alpha_(alpha), zeroIdx_(0) {
				while ((zeroIdx_ < relativeStrikes_.size() - 1) && (relativeStrikes_[zeroIdx_] < 0.0)) ++zeroIdx_;
				QL_REQUIRE(zeroIdx_ > 0, "zeroIdx_ > 0 required");
			}

			virtual Real value(const Array& x) const {
			    Disposable<Array> y = values(x);
			    Real res = 0.0;
			    for (size_t k = 0; k<y.size(); ++k) res += y[k] * y[k];
			    return res / 2.0;
		    }

			boost::shared_ptr<VanillaLocalVolModel> model(const Array& x) const {
				QL_REQUIRE(x.size() == relativeStrikes_.size() - 1, "(x.size() == relativeStrikes_.size() - 1 required");
				std::vector<Real> Sm(zeroIdx_), Mm(zeroIdx_);
				std::vector<Real> Sp(relativeStrikes_.size() - zeroIdx_), Mp(relativeStrikes_.size() - zeroIdx_);  // we need an additional segment for high-strike extrapolation control (CMS)
				for (size_t k = 0; k < Sm.size(); ++k) {
					Sm[k] = model_->forward() + relativeStrikes_[zeroIdx_ - 1 - k];
					Mm[k] = TemplateAuxilliaries::direct<Real>(x[k], minSlope_, maxSlope_);
					if (k > 0) Mm[k] += Mm[k - 1];  // we want to model relative offsets in slope
				}
				for (size_t k = 0; k < Sp.size() - 1; ++k) {  // the last segment is for high-strike extrapolation
					Sp[k] = model_->forward() + relativeStrikes_[zeroIdx_ + 1 + k];
					Mp[k] = TemplateAuxilliaries::direct<Real>(x[Sm.size() + k], minSlope_, maxSlope_);
					if (k > 0) Mp[k] += Mp[k - 1];  // we want to model relative offsets in slope
				}
				Sp[Sp.size() - 1] = model_->forward() + extrapolationRelativeStrike_;
				Mp[Mp.size() - 1] = extrapolationSlope_;
				// now we can create a new model...
				return boost::shared_ptr<VanillaLocalVolModel>(
					new VanillaLocalVolModel(model_->timeToExpiry(), model_->forward(), model_->sigmaATM(), Sp, Sm, Mp, Mm, model_->maxCalibrationIters(),
						model_->onlyForwardCalibrationIters(), model_->adjustATMFlag(), model_->enableLogging(), model_->useInitialMu(), model_->initialMu()));
			}

			virtual Disposable<Array> values(const Array& x) const {
				boost::shared_ptr<VanillaLocalVolModel> newModel = model(x);
				size_t sizeF = relativeStrikes_.size() - 1; // ATM is already calibrated on-the-fly
				if (alpha_ > 0.0) sizeF += relativeStrikes_.size() - 3;  // we may want to minimize curvature, but exclude skew/convexity around ATM
				                                                         // this relies on having at least three relative strikes provided
				Array objectiveF(sizeF);
				for (size_t k = 0; k < zeroIdx_; ++k) {
					objectiveF[k] = newModel->expectation(false, newModel->forward() + relativeStrikes_[k]);
					objectiveF[k] = (objectiveF[k] - smilePrices_[k]) / vanillaVega_[k];
				}
				for (size_t k = zeroIdx_ + 1; k < relativeStrikes_.size(); ++k) {
					objectiveF[k-1] = newModel->expectation(true, newModel->forward() + relativeStrikes_[k]);
					objectiveF[k-1] = (objectiveF[k-1] - smilePrices_[k]) / vanillaVega_[k];
				}
				if (alpha_ > 0.0) {
					for (size_t k = 0; k < zeroIdx_ - 1; ++k) objectiveF[relativeStrikes_.size() - 1 + k] = alpha_ * TemplateAuxilliaries::direct<Real>(x[k], minSlope_, maxSlope_);
					for (size_t k = zeroIdx_ + 2; k < relativeStrikes_.size(); ++k) objectiveF[(relativeStrikes_.size() - 1) + (zeroIdx_ - 1) + k] = alpha_ * TemplateAuxilliaries::direct<Real>(x[k-1], minSlope_, maxSlope_);
				}
				return objectiveF;
			}

			Disposable<Array> initialValues() const { // we use zero slope as initial guess
				return Array(relativeStrikes_.size(), TemplateAuxilliaries::inverse<Real>(0.0, minSlope_, maxSlope_));
			}

		};

	}

	VanillaLocalVolModelSmileSection::VanillaLocalVolModelSmileSection(
		const Date&                                       expiryDate,
		const Rate&                                       forward,
		const std::vector<Rate>&                          relativeStrikes,
		const std::vector<Volatility>&                    smileVolatilities,
		const Real                                        extrapolationRelativeStrike,
		const Real                                        extrapolationSlope,
		bool                                              vegaWeighted,
		const boost::shared_ptr<EndCriteria>&             endCriteria,
		const boost::shared_ptr<OptimizationMethod>&      method,
		const DayCounter&                                 dc,
		const Date&                                       referenceDate,
		const VolatilityType                              type,
		const Rate                                        shift,
		const boost::shared_ptr<VanillaLocalVolModel>&    model,
		const Real                                        minSlope,    //  lower boundary for m in calibration
		const Real                                        maxSlope,    //  upper boundary for m in calibration
		const Real                                        alpha )      //  Tikhonov alpha
		: SmileSection(expiryDate, dc, referenceDate, type, shift) {
		// first we check for consistent strike inputs...
		QL_REQUIRE(relativeStrikes.size() >= 3, "relativeStrikes.size() >= 3 required");
		for (size_t k = 1; k < relativeStrikes.size(); ++k) QL_REQUIRE(relativeStrikes[k] > relativeStrikes[k - 1], "relativeStrikes[k] > relativeStrikes[k - 1] required");
		size_t zeroIdx = 0;
		while ((zeroIdx < relativeStrikes.size() - 1) && (relativeStrikes[zeroIdx] < 0.0)) ++zeroIdx;
		QL_REQUIRE(relativeStrikes[zeroIdx] == 0.0, "relativeStrikes[zeroIdx] == 0.0 required");
		QL_REQUIRE(zeroIdx > 0, "zeroIdx > 0 required");
		QL_REQUIRE(zeroIdx < relativeStrikes.size() - 1, "zeroIdx < relativeStrikes.size() - 1 required");
		// now we can check for consistent volatility inputs
		QL_REQUIRE(smileVolatilities.size() == relativeStrikes.size(), "smileVolatilities.size() == relativeStrikes.size() required");
		for (size_t k = 0; k < smileVolatilities.size(); ++k) QL_REQUIRE(smileVolatilities[k] > 0.0, "smileVolatilities[k] > 0.0 required");
		// calibration should be independent of the volatility type, thus we calculate prices and vega based on input vol-type
		Real timeToExpiry = dc.yearFraction(this->referenceDate(), expiryDate);
		QL_REQUIRE(timeToExpiry > 0.0, "timeToExpiry > 0.0 required");
		std::vector<Real> smilePrices(relativeStrikes.size()), vanillaVega(relativeStrikes.size());
		for (size_t k = 0; k < relativeStrikes.size(); ++k) {
			if (type == Normal) {
				smilePrices[k] = bachelierBlackFormula((relativeStrikes[k] < 0.0) ? (Option::Put) : (Option::Call), forward + relativeStrikes[k], forward, smileVolatilities[k] * sqrt(timeToExpiry));
				if (vegaWeighted) vanillaVega[k] = 1.0;  // since we already calibrate to prices we don't want additional Vega weighting
				else vanillaVega[k] = bachelierBlackFormulaStdDevDerivative(forward + relativeStrikes[k], forward, smileVolatilities[k] * sqrt(timeToExpiry)) * sqrt(timeToExpiry);
			}
			else {  // ShiftedLognormal
				smilePrices[k] = blackFormula((relativeStrikes[k] < 0.0) ? (Option::Put) : (Option::Call), forward + relativeStrikes[k], forward, smileVolatilities[k] * sqrt(timeToExpiry), 1.0, shift);
				if (vegaWeighted) vanillaVega[k] = 1.0;  // since we already calibrate to prices we don't want additional Vega weighting
				else vanillaVega[k] = blackFormulaStdDevDerivative(forward + relativeStrikes[k], forward, smileVolatilities[k] * sqrt(timeToExpiry),1.0,shift) * sqrt(timeToExpiry);
			}
		}
		// we also need a consistent model for calibration, thus we set up a normal model matching ATM prices
		Real sigmaATM = bachelierBlackFormulaImpliedVol(Option::Call, forward, forward, timeToExpiry, smilePrices[zeroIdx]);
		std::vector<Real> Sm(1, forward + relativeStrikes[zeroIdx - 1]), Sp(1, forward + relativeStrikes[zeroIdx - 1]), Mm(1, 0.0), Mp(1, 0.0);
		if (model) model_ = boost::shared_ptr<VanillaLocalVolModel>(
			new VanillaLocalVolModel(timeToExpiry, forward, sigmaATM, Sp, Sm, Mp, Mm, model->maxCalibrationIters(), model->onlyForwardCalibrationIters(), model->adjustATMFlag(), model->enableLogging(), model->useInitialMu(), model->initialMu()));
		else model_ = boost::shared_ptr<VanillaLocalVolModel>(
			new VanillaLocalVolModel(timeToExpiry, forward, sigmaATM, Sp, Sm, Mp, Mm, 100, 0, true, true, false, 0.0));
		// now we may set up the optimization problem...
		VanillaLocalVolHelper costFunction(model_, relativeStrikes, smilePrices, vanillaVega, extrapolationRelativeStrike, extrapolationSlope, minSlope, maxSlope, alpha);
		Problem problem(costFunction, NoConstraint(), costFunction.initialValues());
		method->minimize(problem, *endCriteria);
		model_ = costFunction.model(problem.currentValue());
	}



	Real VanillaLocalVolModelSmileSection::volatilityImpl(Rate strike) const {
		Option::Type type = (strike >= model_->forward()) ? (Option::Call) : (Option::Put);
		Real price = model_->expectation((type == Option::Call) ? (true) : (false), strike);
		if (volatilityType() == Normal) return bachelierBlackFormulaImpliedVol(type, strike, model_->forward(), model_->timeToExpiry(), price);
		else return blackFormulaImpliedStdDev(type, strike, model_->forward(), price);
	}

	Real VanillaLocalVolModelSmileSection::optionPrice(Rate strike, Option::Type type, Real discount) const {
		Real otmPrice = model_->expectation((strike>= model_->forward()) ? (true) : (false), strike);
		if ((type == Option::Call) && (strike < model_->forward())) return (otmPrice + (model_->forward() - strike))*discount;
		if ((type == Option::Put) && (strike > model_->forward())) return (otmPrice - (model_->forward() - strike))*discount;
		return otmPrice*discount;
	}


}

