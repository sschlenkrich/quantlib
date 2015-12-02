/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_template2fnormalmodel_hpp
#define quantlib_template2fnormalmodel_hpp

#include <boost/shared_ptr.hpp>
#include <ql/errors.hpp>
//#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>
//#include <ql/experimental/template/auxilliaries/solver1d.hpp>
#include <ql/experimental/template/templatestochasticprocess.hpp>



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

		Handle<AssetTermStructure> phi_;            // deterministic part

		// term structure for deterministic part
		
		// unique grid for time-dependent parameters
		VecD                       times_;   // time-grid of left-constant model parameter values
		// time-dependent parameters, left-piecewise constant on times_-grid
		VecA                       sigma_;  // volatility for Y
		VecA                       eta_;    // volatility for Z
		// time-homogeneous parameters
		PassiveType                a_;      // mean reversion for Y 		
		PassiveType                b_;      // mean reversion for Z
		PassiveType                rho_;    // correlation Y vs Z
		
	public:
		// constructor
		Template2FNormalModel( const Handle<AssetTermStructure> &    phi,
		                       const VecD &                          times,
							   const VacA &                          sigma,
							   const VecA &                          eta,
							   const PassiveType                     a,
							   const PassiveType                     b,
							   const PassiveType                     rho )
		: phi_(phi), times_(times), sigma_(sigma), eta_(eta), a_(a), b_(b),	rho_(rho) {
			// check for valid parameter inputs
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
			// TODO...
			// Y-variable z(t) (S(t) + lambda)^beta dW(t)
			B[0][0] = X[1] * pow(X[0] + lambda_, beta_);
			B[0][1] = 0.0;
			// z-variable nu z(t) dZ(t)
			B[1][0] = rho_ * nu_*X[1];
			B[1][1] = sqrt(1-rho_*rho_) * nu_*X[1];
			// finished
			return B;
		}


	};

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2fnormalmodel_hpp */
