/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_template2fnormalmodel_hpp
#define quantlib_template2fnormalmodel_hpp

#include <boost/shared_ptr.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>
#include <ql/experimental/template/auxilliaries/templateintegrators.hpp>
//#include <ql/experimental/template/auxilliaries/solver1d.hpp>
#include <ql/experimental/template/templatestochasticprocess.hpp>
#include <ql/experimental/template/commodity/indextermstructure.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // 2-factor mean reverting model
	//
	//     X(t) = phi(t) + Y(t) + Z(t), X(0)=Y(0)=0
	//    dX(t) = -a X(t) dt  +  sigma(t) dW_Y(t)
	//    dY(t) = -b Y(t) dt  +    eta(t) dW_Z(t)
	//    dW_Y(t) dW_Z(t) = rho dt
	//
	template <class DateType, class PassiveType, class ActiveType>
	class Template2FNormalModel : public TemplateStochasticProcess<DateType,PassiveType,ActiveType> {
	protected:
	
		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
	
		// attributes defining the model
		Handle<YieldTermStructure> termStructure_;  // NOT USED; the yield curve is assumed to be passive

		// term structure for deterministic part
		Handle<IndexTermStructure> phi_;            // deterministic part

		
		// unique grid for time-dependent parameters
		VecD                       times_;   // time-grid of left-constant model parameter values
		// time-dependent parameters, left-piecewise constant on times_-grid
		VecA                       sigma_;  // volatility for Y
		VecA                       eta_;    // volatility for Z
		// time-homogeneous parameters
		PassiveType                a_;      // mean reversion for Y 		
		PassiveType                b_;      // mean reversion for Z
		PassiveType                rho_;    // correlation Y vs Z

		inline const ActiveType sigma(const DateType t) const { return sigma_[TemplateAuxilliaries::idx(times_,t)]; }
		inline const ActiveType eta  (const DateType t) const { return eta_[TemplateAuxilliaries::idx(times_,t)];   }

		// functor for (co-)variance integration
		struct CovarianceHelper {
			ActiveType a_, b_;
			DateType T_;
			CovarianceHelper(const ActiveType a, const ActiveType b, const DateType T) : a_(a), b_(b), T_(T) {}
			// f(t) = exp{-(a+b)(T-t)}
			// F(t) = exp{-(a+b)(T-t)} / (a+b)
			inline PassiveType operator() (const DateType t) { return exp(-(a_+b_)*(T_-t)) / (a_+b_); }
		};

	public:
		// constructor
		
		Template2FNormalModel( const Handle<IndexTermStructure>&    phi,
		                       const VecD&                          times,
							   const VecA&                          sigma,
							   const VecA&                          eta,
							   const PassiveType                     a,
							   const PassiveType                     b,
							   const PassiveType                     rho 
							   )
		: phi_(phi), times_(times), sigma_(sigma), eta_(eta), a_(a), b_(b),	rho_(rho) {
			// check for valid parameter inputs
		}
	
		// analytic formulas

		// future expectation
        inline ActiveType futureAsset(const DateType t, const DateType T, const ActiveType Y, const ActiveType Z) {
			return phi_->value(T) + exp(-a_*(T-t))*Y + exp(-b_*(T-t))*Z;
		}

		// (future) variance of Y process
		inline ActiveType varianceY(const DateType t, const DateType T) {
			CovarianceHelper F(a_,a_,T);
			VecA sigma2(sigma_.size());
			for (size_t k=0; k<sigma_.size(); ++k) sigma2[k] = sigma_[k]*sigma_[k];
			TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, sigma2, F);
			ActiveType variance = integral(t,T);
			return variance;
		}

		// (future) variance of Z process
		inline ActiveType varianceZ(const DateType t, const DateType T) {
			CovarianceHelper F(b_,b_,T);
			VecA eta2(eta_.size());
			for (size_t k=0; k<sigma_.size(); ++k) eta2[k] = eta_[k]*eta_[k];
			TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, eta2, F);
			ActiveType variance = integral(t,T);
			return variance;
		}

		// (future) covariance X-Z process
		inline ActiveType covarianceYZ(const DateType t, const DateType T) {
			CovarianceHelper F(a_,b_,T);
			VecA sigmaTimesEta(sigma_.size());
			for (size_t k=0; k<sigmaTimesEta.size(); ++k) sigmaTimesEta[k] = sigma_[k]*eta_[k];
			TemplateAuxilliaries::PieceWiseConstantIntegral<PassiveType,ActiveType,CovarianceHelper> integral(times_, sigmaTimesEta, F);
			ActiveType covariance = integral(t,T);
			return covariance;
		}

		// variance of asset future
		inline ActiveType varianceFuture( const DateType startTime, const DateType expiryTime, const DateType settlementTime ) {
			PassiveType B = exp(-a_*(settlementTime-expiryTime));
			PassiveType C = exp(-b_*(settlementTime-expiryTime));
			ActiveType var = B*B*varianceY(startTime,expiryTime) + C*C*varianceZ(startTime,expiryTime) + 2.0*rho_*B*C*covarianceYZ(startTime,expiryTime);
			return var;
		}

		// basic instruments

        inline ActiveType averageFuture ( const VecD& settlementTimes, const VecP& settlementWeights) {
			ActiveType fut = 0.0;
			size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
			for (size_t k=0; k<N; ++k) fut += settlementWeights[k]*phi_->value(settlementTimes[k]);
			return fut;
		}

        inline ActiveType varianceAverageFuture ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights) {
			PassiveType B = 0, C = 0;
			size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
			for (size_t k=0; k<N; ++k) {
				B += settlementWeights[k]*exp(-a_*(settlementTimes[k]-expiryTime));
				C += settlementWeights[k]*exp(-b_*(settlementTimes[k]-expiryTime));
			}
			ActiveType var = B*B*varianceY(0.0,expiryTime) + C*C*varianceZ(0.0,expiryTime) + 2.0*rho_*B*C*covarianceYZ(0.0,expiryTime);
			return var;
		}

		inline ActiveType vanillaOption ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, PassiveType strike, int callOrPut) {
			ActiveType fut = averageFuture(settlementTimes,settlementWeights);
			ActiveType var = varianceAverageFuture(expiryTime,settlementTimes,settlementWeights);
			ActiveType pv  = TemplateAuxilliaries::Bachelier(fut,strike,sqrt(var),1.0,callOrPut);
			return pv;
		}

		// stochastic process interface

		// dimension of X
		inline virtual size_t size() { return 2; }
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline virtual size_t factors() { return 2; }
		// initial values for simulation
		inline virtual VecP initialValues() {
			VecP X(2,0.0);
			return X;
		}
		// a[t,X(t)]
		inline virtual VecA drift( const DateType t, const VecA& X) {
			VecA a(2);
			// Y-variable
			a[0] = -a_ * X[0];
            // Z-variable
			a[1] = -b_ * X[1];
			return a;
		}
		// b[t,X(t)]
		inline virtual MatA diffusion( const DateType t, const VecA& X) {
			MatA B(2);
			B[0].resize(2);
			B[1].resize(2);
			// Y-variable sigma(t) dW_Y(t)
			B[0][0] = sigma(t); 
			B[0][1] = 0.0;
			// Z-variable eta(t) dW_Z(t)
			B[1][0] = eta(t) * rho_             ;
			B[1][1] = eta(t) * sqrt(1-rho_*rho_);
			// finished
			return B;
		}

		// stochastic process variables and payoffs

		inline virtual ActiveType asset(const DateType t, const VecA& X)                    { return futureAsset(t,t,X[0],X[1]);  }
		inline virtual ActiveType future(const DateType t, const DateType T, const VecA& X) { return futureAsset(t,T,X[0],X[1]);  }

	};

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2fnormalmodel_hpp */
