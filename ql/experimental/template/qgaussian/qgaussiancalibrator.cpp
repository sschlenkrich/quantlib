/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/


#include <ql/experimental/template/qgaussian/qgaussiancalibrator.hpp>

namespace QuantLib {

	QuasiGaussianModelCalibrator::Objective::Objective(
		          QuasiGaussianModelCalibrator              *calibrator,
			      const std::vector< std::vector< bool > >&  isInput,
				  std::vector< std::vector< bool > >&        isOutput ) 
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
		swaptions_.resize(outputSize_);
		size_t idx = 0;
		for (size_t i=0; i<isOutput_.size(); ++i) {
			for (size_t j=0; j<isOutput_[i].size(); ++j) {
				if (isOutput_[i][j]) {
					swaptions_[idx] = boost::shared_ptr<CalibSwaption>(new CalibSwaption(calibrator_->lambda_[i][j],calibrator_->b_[i][j],calibrator_->eta_[i][j],calibrator_->swaptions_[i][j],model_->termStructure(), true) );
					++idx;
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

	void QuasiGaussianModelCalibrator::Objective::update(const Array& X) {
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

}

