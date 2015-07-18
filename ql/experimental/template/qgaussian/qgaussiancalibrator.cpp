/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/


#include <ql/experimental/template/qgaussian/qgaussiancalibrator.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>


namespace QuantLib {

	QuasiGaussianModelCalibrator::Objective::Objective(
		          QuasiGaussianModelCalibrator              *calibrator,
			      const std::vector< std::vector< bool > >&  isInput,
				  const std::vector< std::vector< bool > >&  isOutput ) 
	    : calibrator_(calibrator), isInput_(isInput), isOutput_(isOutput) {
		// copy model initial values
		model_ = boost::shared_ptr<RealQuasiGaussianModel>(new RealQuasiGaussianModel(*calibrator->model_));
		// checking dimensions
	    QL_REQUIRE(isInput_.size()==model_->times().size(),"QuasiGaussianModelCalibrator::Objective: wrong input dimension.");
		for (size_t i=0; i<isInput_.size(); ++i) {
			QL_REQUIRE(isInput_[i].size()==(2*model_->factors()-1),"QuasiGaussianModelCalibrator::Objective: wrong input dimension.");
		}
		QL_REQUIRE(isOutput_.size()==calibrator_->swaptions_.size(),"QuasiGaussianModelCalibrator::Objective: wrong output dimension.");
		for (size_t i=0; i<isOutput_.size(); ++i) {
			QL_REQUIRE(isOutput_[i].size()==(3*calibrator_->swaptions_[i].size()),"QuasiGaussianModelCalibrator::Objective: wrong output dimension.");
		}
		// count inputs and outputs
		inputSize_=0;
		for (size_t i=0; i<isInput_.size(); ++i)
			for (size_t j=0; j<isInput_[i].size(); ++j) if (isInput_[i][j]) ++inputSize_;
		outputSize_=0;
		for (size_t i=0; i<isOutput_.size(); ++i)
			for (size_t j=0; j<isOutput_[i].size(); ++j) if (isOutput_[i][j]) ++outputSize_;

		// set up target swaptions
		// we may have a calibration swaption for each input swaption
		calibSwaptions_.resize(calibrator_->swaptions_.size());
		for (size_t i=0; i<calibSwaptions_.size(); ++i) {
			calibSwaptions_[i].resize(calibrator_->swaptions_[i].size());
			for (size_t j=0; j<calibSwaptions_[i].size(); ++j) {
				// allocate only if at least one of lambda, b or eta is calibrated
				if (isOutput_[i][j]                                   ||
					isOutput_[i][j+calibrator_->swaptions_[i].size()] || 
					isOutput_[i][j+2*calibrator_->swaptions_[i].size()] ) {
					// assume equal dimension of lambda, b, eta and swaptions
					calibSwaptions_[i][j] = boost::shared_ptr<CalibSwaption>(new CalibSwaption(calibrator_->lambda_[i][j],calibrator_->b_[i][j],calibrator_->eta_[i][j],calibrator_->swaptions_[i][j],model_->termStructure(), true) );
				}
			}
		}
		// ready for optimisation...
	}

	Array QuasiGaussianModelCalibrator::Objective::initialise() {
		Array X(inputSize_);
		size_t d = model_->factors()-1;
		size_t idx = 0;
		for (size_t i=0; i<isInput_.size(); ++i) {
			for (size_t j=0; j<d; ++j) {
				if (isInput_[i][j]) {
					Real lambda = model_->lambda()[j][i];
					Real x = QuasiGaussianModelCalibrator::inverse(lambda,calibrator_->lambdaMin_,calibrator_->lambdaMax_);
					X[idx] = x;
					++idx;
				}
			}
			for (size_t j=0; j<d; ++j) {
				if (isInput_[i][d+j]) {
					Real b = model_->b()[j][i];
					Real x = QuasiGaussianModelCalibrator::inverse(b,calibrator_->bMin_,calibrator_->bMax_);
					X[idx] = x;
					++idx;
				}
			}
			if (isInput_[i][d+d]) {
				Real eta = model_->eta()[i];
				Real x = QuasiGaussianModelCalibrator::inverse(eta,calibrator_->etaMin_,calibrator_->etaMax_);
				X[idx] = x;
				++idx;
			}
		}
		return X;
	}

	void QuasiGaussianModelCalibrator::Objective::update(const Array& X) const {
		std::vector< std::vector< Real > > m_lambda = model_->lambda();
		std::vector< std::vector< Real > > m_b      = model_->b();
		std::vector< Real >                m_eta    = model_->eta();
		size_t d = model_->factors()-1;
		size_t idx = 0;
		for (size_t i=0; i<isInput_.size(); ++i) {
			for (size_t j=0; j<d; ++j) {
				if (isInput_[i][j]) {
					Real x = X[idx];
					Real lambda = QuasiGaussianModelCalibrator::direct(x,calibrator_->lambdaMin_,calibrator_->lambdaMax_);
					m_lambda[j][i] = lambda;
					++idx;
				}
			}
			for (size_t j=0; j<d; ++j) {
				if (isInput_[i][d+j]) {
					Real x = X[idx];
					Real b = QuasiGaussianModelCalibrator::direct(x,calibrator_->bMin_,calibrator_->bMax_);
					m_b[j][i] = b;
					++idx;
				}
			}
			if (isInput_[i][d+d]) {
				Real x = X[idx];
				Real eta = QuasiGaussianModelCalibrator::direct(x,calibrator_->etaMin_,calibrator_->etaMax_);
				m_eta[i] = eta;
				++idx;
			}
		}
		model_->update(m_lambda,m_b,m_eta);
	}

	Disposable<Array> QuasiGaussianModelCalibrator::Objective::values(const Array& x) const {
		update(x);
		// we may have a swaption model for each swaption
		std::vector< std::vector< boost::shared_ptr<RealQGSwaptionModel> > > swaptionModels;
		// however we build the model only if neccessary
		swaptionModels.resize(calibrator_->swaptions_.size());
		for (size_t i=0; i<swaptionModels.size(); ++i) {
			swaptionModels[i].resize(calibrator_->swaptions_[i].size());
			for (size_t j=0; j<swaptionModels[i].size(); ++j) {
				if (calibSwaptions_[i][j]) {										
					swaptionModels[i][j] = boost::shared_ptr<RealQGSwaptionModel>( new RealQGSwaptionModel( model_, calibSwaptions_[i][j]->floatTimes(), calibSwaptions_[i][j]->floatWeights(), calibSwaptions_[i][j]->fixedTimes(), calibSwaptions_[i][j]->fixedWeights(),calibrator_->modelTimes_, calibrator_->useExpectedXY_ ) );
				}
			}
		}
		// averaging and objective calculation
		Array objective(outputSize_);
		size_t idx=0;
		for (size_t i=0; i<swaptionModels.size(); ++i) {
			for (size_t j=0; j<swaptionModels[i].size(); ++j) {
				if (isOutput_[i][j]) {
					objective[idx] = swaptionModels[i][j]->averageLambda( calibSwaptions_[i][j]->exerciseTimes()[0] );
					objective[idx] -= calibSwaptions_[i][j]->lambda();
					++idx;
				}
				if (isOutput_[i][j+calibSwaptions_[i].size()]) {
					objective[idx] = swaptionModels[i][j]->averageB( calibSwaptions_[i][j]->exerciseTimes()[0] );
					objective[idx] -= calibSwaptions_[i][j]->b();
					++idx;
				}
				if (isOutput_[i][j+2*calibSwaptions_[i].size()]) {
					objective[idx] = swaptionModels[i][j]->averageEta( calibSwaptions_[i][j]->exerciseTimes()[0] );
					objective[idx] -= calibSwaptions_[i][j]->eta();
					++idx;
				}
			}
		}
		return objective;
	}

	Real QuasiGaussianModelCalibrator::Objective::value(const Array& x) const {
		Real sum = 0.0;
		for (size_t i=0; i<x.size(); ++i) sum += x[i]*x[i];
		return sum;
	}

	const boost::shared_ptr<RealQuasiGaussianModel> QuasiGaussianModelCalibrator::Objective::model(const Array& x) {
		update(x);
		return model_;
	}

	// constructor
	QuasiGaussianModelCalibrator::QuasiGaussianModelCalibrator(
		                              boost::shared_ptr<RealQuasiGaussianModel> model,
			                          boost::shared_ptr<RealMCSimulation>       mcSimulation,
									  std::vector< std::vector< boost::shared_ptr<Swaption> > > swaptions,
									  std::vector< std::vector< Real > > lambda,
                                      std::vector< std::vector< Real > > b,
                                      std::vector< std::vector< Real > > eta,
									  Real                               lambdaMin,
									  Real                               lambdaMax,
									  Real                               bMin,
									  Real                               bMax,
									  Real                               etaMin,
									  Real                               etaMax,
									  std::vector< Real >                modelTimes,
                                      bool                               useExpectedXY) 
									  : model_(model), mcSimulation_(mcSimulation), swaptions_(swaptions), lambda_(lambda), b_(b), eta_(eta),
									    lambdaMin_(lambdaMin), lambdaMax_(lambdaMax), bMin_(bMin), bMax_(bMax), etaMin_(etaMin), etaMax_(etaMax),
										modelTimes_(modelTimes), useExpectedXY_(useExpectedXY) {
        // check dimensions
		QL_REQUIRE(swaptions_.size()>0,"QuasiGaussianModelCalibrator: Wrong swaptions dimension." );
		for (size_t i=0; i<swaptions_.size(); ++i) QL_REQUIRE(swaptions_[i].size()>0,"QuasiGaussianModelCalibrator: Wrong swaptions dimension." );

		QL_REQUIRE(lambda_.size()==swaptions_.size(),"QuasiGaussianModelCalibrator: Wrong lambda dimension." );
		for (size_t i=0; i<lambda_.size(); ++i) QL_REQUIRE(lambda_[i].size()==swaptions_[i].size(),"QuasiGaussianModelCalibrator: Wrong lambda dimension." );

		QL_REQUIRE(b_.size()==swaptions_.size(),"QuasiGaussianModelCalibrator: Wrong b dimension." );
		for (size_t i=0; i<b_.size(); ++i) QL_REQUIRE(b_[i].size()==swaptions_[i].size(),"QuasiGaussianModelCalibrator: Wrong b dimension." );

		QL_REQUIRE(eta_.size()==swaptions_.size(),"QuasiGaussianModelCalibrator: Wrong eta dimension." );
		for (size_t i=0; i<eta_.size(); ++i) QL_REQUIRE(eta_[i].size()==swaptions_[i].size(),"QuasiGaussianModelCalibrator: Wrong eta dimension." );

		// better check and adjust model times...
	}

	const boost::shared_ptr<RealQuasiGaussianModel> QuasiGaussianModelCalibrator::calibrate(
						           const std::vector< std::vector< bool > >&  isInput,
			                       const std::vector< std::vector< bool > >&  isOutput,
								   	// optimization parameters
		                           Real                                       epsfcn,
								   Real                                       ftol,  
								   Real                                       xtol,  
								   Real                                       gtol,  
		                           Size                                       maxfev ) {
		// set up
		NoConstraint constraint;
		Objective obj(this, isInput, isOutput);
		Array x = obj.initialise();
		Problem problem(obj, constraint, x);
		LevenbergMarquardt optimizationMethod(epsfcn, xtol, gtol);
		EndCriteria endCriteria(maxfev, 100 /* unused */, 0 /* unused */, ftol, 0 /* unused */);
		// calibrate
		optimizationMethod.minimize(problem,endCriteria);
		return obj.model(problem.currentValue());
	    return model_;

	}


}

