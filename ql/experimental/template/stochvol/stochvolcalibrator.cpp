/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

#include <ql/pricingengines/blackformula.hpp>

#include <ql/experimental/template/stochvol/stochvolcalibrator.hpp>

namespace QuantLib {

		// count number of effective states
		size_t StochVolModelCalibrator::xDim() const {
			size_t N=0;
			if (!lambdaIsFixed_ ) ++N;
			if (!bIsFixed_		) ++N;
			if (!LIsFixed_		) ++N;
			if (!thetaIsFixed_	) ++N;
			if (!mIsFixed_		) ++N;
			if (!etaIsFixed_	) ++N;
			if (!z0IsFixed_		) ++N;
			if (!rhoIsFixed_	) ++N;
			return N;
		}

		Array StochVolModelCalibrator::direct(const Array& X) const {
			QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
			Array Y(X.size());
			size_t idx = 0;
			if (!lambdaIsFixed_ ) { Y[idx] = 2*(atan(X[idx])/M_PI + 0.5);  ++idx; }   // 0 < lambda < 2
			if (!bIsFixed_		) { Y[idx] = atan(X[idx])/M_PI + 0.5;      ++idx; }   // 0 < b < 1
			if (!LIsFixed_		) { Y[idx] = exp(X[idx]);                  ++idx; }   // L > 0
			if (!thetaIsFixed_	) { Y[idx] = exp(X[idx]);                  ++idx; }   // theta > 0
			if (!mIsFixed_		) { Y[idx] = exp(X[idx]);                  ++idx; }   // m > 0 
			if (!etaIsFixed_	) { Y[idx] = 1*(atan(X[idx])/M_PI + 0.5);  ++idx; }   // 0 < eta < 1
			if (!z0IsFixed_		) { Y[idx] = exp(X[idx]);                  ++idx; }   // z0 > 0
			if (!rhoIsFixed_	) { Y[idx] = atan(X[idx])*2.0/M_PI;        ++idx; }   // -1 < rho < 1
			return Y;
		}

		Array StochVolModelCalibrator::inverse(const Array& Y) const {
			QL_REQUIRE(Y.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
			Array X(Y.size());
			size_t idx = 0;
			if (!lambdaIsFixed_ ) { X[idx] = tan((Y[idx]/2-0.5)*M_PI);  ++idx; }
			if (!bIsFixed_		) { X[idx] = tan((Y[idx]-0.5)*M_PI);    ++idx; } 
			if (!LIsFixed_		) { X[idx] = log(Y[idx]);               ++idx; } 
			if (!thetaIsFixed_	) { X[idx] = log(Y[idx]);               ++idx; } 
			if (!mIsFixed_		) { X[idx] = log(Y[idx]);               ++idx; } 
			if (!etaIsFixed_	) { X[idx] = tan((Y[idx]/1-0.5)*M_PI);  ++idx; }
			if (!z0IsFixed_		) { X[idx] = log(Y[idx]);               ++idx; }
			if (!rhoIsFixed_	) { X[idx] = tan(Y[idx]*M_PI/2.0);      ++idx; }
			return X;
		}

		// initialize state X with model parameters and apply inverse transformation
		Array StochVolModelCalibrator::initialise() {
			Array X(xDim());
			size_t idx = 0;
			if (!lambdaIsFixed_ ) { X[idx] = lambda_;   ++idx; }
			if (!bIsFixed_		) { X[idx] = b_;		++idx; } 
			if (!LIsFixed_		) { X[idx] = L_;		++idx; } 
			if (!thetaIsFixed_	) { X[idx] = theta_;	++idx; } 
			if (!mIsFixed_		) { X[idx] = m_;		++idx; } 
			if (!etaIsFixed_	) { X[idx] = eta_;	    ++idx; }
			if (!z0IsFixed_		) { X[idx] = z0_;	    ++idx; }
			if (!rhoIsFixed_	) { X[idx] = rho_;	    ++idx; }
			return inverse(X);
		}

		// create model for given state vector
		boost::shared_ptr<RealStochVolModel> StochVolModelCalibrator::model(const Array& X) const {
		    // model parameters
		    Real  lambda;
		    Real  b;
		    Real  L;
		    Real  theta;
		    Real  m;
		    Real  eta;
		    Real  z0;
		    Real  rho;
			// get parameters
			QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
			size_t idx = 0;
			if (!lambdaIsFixed_ ) { lambda = X[idx];  ++idx; } else { lambda = lambda_; }
			if (!bIsFixed_		) { b      = X[idx];  ++idx; } else { b      = b_;     	}
			if (!LIsFixed_		) { L      = X[idx];  ++idx; } else { L      = L_;     	}
			if (!thetaIsFixed_	) { theta  = X[idx];  ++idx; } else { theta  = theta_; 	}
			if (!mIsFixed_		) { m	   = X[idx];  ++idx; } else { m	     = m_;	   	}
			if (!etaIsFixed_	) { eta    = X[idx];  ++idx; } else { eta    = eta_;    }
			if (!z0IsFixed_		) { z0     = X[idx];  ++idx; } else { z0     = z0_;     }
			if (!rhoIsFixed_	) { rho    = X[idx];  ++idx; } else { rho    = rho_;    }
			return boost::shared_ptr<RealStochVolModel>(new RealStochVolModel(lambda,b,L,theta,m,eta,z0,rho));
		}

		Disposable<Array> StochVolModelCalibrator::values(const Array& x) const {
			boost::shared_ptr<RealStochVolModel> m = model(direct(x));
			std::vector<Real>          fwPrices(strikes_.size());
			std::vector<Option::Type>  callOrPut(strikes_.size());
			std::vector<Real>          normalVols(strikes_.size());
			Array                      volVariance(strikes_.size());
			for (size_t k=0; k<strikes_.size(); ++k) {
				callOrPut[k]   = (strikes_[k] < forward_) ? (Option::Put) : (Option::Call);
				fwPrices[k]    = m->vanillaOption(forward_,strikes_[k],exercTime_,callOrPut[k],1.0e-12,10000);
				normalVols[k]  = bachelierBlackFormulaImpliedVol(callOrPut[k],strikes_[k],forward_,exercTime_,fwPrices[k]);
				volVariance[k] = normalVols[k] - vols_[k];
			}
			return volVariance;
		}

        Real StochVolModelCalibrator::value(const Array& x) const {
			Disposable<Array> y = values(x);
			Real res=0.0;
			for (size_t k=0; k<y.size(); ++k) res += y[k]*y[k];
			return res/2.0;
		}

		// update model parameters with state X
		void StochVolModelCalibrator::update(const Array& X) {
			QL_REQUIRE(X.size()==xDim(),"StochVolModelCalibrator Error: dimension mismatch.");
			size_t idx = 0;
			if (!lambdaIsFixed_ ) { lambda_ = X[idx];  ++idx; }
			if (!bIsFixed_		) { b_      = X[idx];  ++idx; } 
			if (!LIsFixed_		) { L_      = X[idx];  ++idx; } 
			if (!thetaIsFixed_	) { theta_  = X[idx];  ++idx; } 
			if (!mIsFixed_		) { m_	    = X[idx];  ++idx; } 
			if (!etaIsFixed_	) { eta_    = X[idx];  ++idx; }
			if (!z0IsFixed_		) { z0_     = X[idx];  ++idx; }
			if (!rhoIsFixed_	) { rho_    = X[idx];  ++idx; }
		}

		// construct object and calibrate
		StochVolModelCalibrator::StochVolModelCalibrator( // initial model parameters
		                         const Real  lambda,
		                         const Real  b,
		                         const Real  L,
		                         const Real  theta,
		                         const Real  m,
		                         const Real  eta,
		                         const Real  z0,
		                         const Real  rho,
		                         // calibration flags
		                         const bool  lambdaIsFixed,
		                         const bool  bIsFixed,
		                         const bool  LIsFixed,
		                         const bool  thetaIsFixed,
		                         const bool  mIsFixed,
		                         const bool  etaIsFixed,
		                         const bool  z0IsFixed,
		                         const bool  rhoIsFixed,
		                         // calibration targets
		                         const Real  exercTime,
		                         const Real  forward,
		                         const std::vector<Real>& strikes,
		                         const std::vector<Real>& vols )
        : lambda_(lambda), b_(b), L_(L), theta_(theta), m_(m), eta_(eta), z0_(z0), rho_(rho),
		  lambdaIsFixed_(lambdaIsFixed), bIsFixed_(bIsFixed), LIsFixed_(LIsFixed), thetaIsFixed_(thetaIsFixed), mIsFixed_(mIsFixed), etaIsFixed_(etaIsFixed), z0IsFixed_(z0IsFixed), rhoIsFixed_(rhoIsFixed),
		  exercTime_(exercTime), forward_(forward), strikes_(strikes), vols_(vols) {
			// set up
			NoConstraint constraint;
			Array x = initialise();
			{
			    Problem problem(*this, constraint, x);
				LevenbergMarquardt optimizationMethod(1e-4, 1e-4, 1e-4);
				EndCriteria endCriteria(60000, 100, 1e-4, 1e-4, 1e-4);
				// calibrate
				optimizationMethod.minimize(problem,endCriteria);
				update(direct(problem.currentValue()));
			}
		}


}

