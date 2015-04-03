/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_templateshiftedsabrmodel_hpp
#define quantlib_templateshiftedsabrmodel_hpp

//#include <complex>
#include <boost/shared_ptr.hpp>
//#include <boost/function.hpp>
#include <ql/errors.hpp>
//#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>
//#include <ql/experimental/template/auxilliaries/gausslobatto.hpp>
//#include <ql/experimental/template/auxilliaries/Complex.hpp>
//#include <ql/experimental/template/auxilliaries/solver1d.hpp>
#include <ql/experimental/template/templatestochasticprocess.hpp>
//#include <ql/experimental/template/stochvol/templatehestonmodel.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // SABR model with shift
	//
	//    dS(t) = z(t) (S(t) + lambda)^beta dW(t)
	//    dz(t) = nu z(t) dZ(t)
	//     z(0) = alpha
	//    dW(t) dZ(t) = rho dt
	//
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateShiftedSABRModel : public TemplateStochasticProcess<DateType,PassiveType,ActiveType> {
	private:
	    ActiveType S0_, lambda_, alpha_, beta_, rho_, nu_;
	public:
		// constructor
		TemplateShiftedSABRModel( ActiveType S0, 
			                      ActiveType lambda, 
								  ActiveType alpha, 
								  ActiveType beta,
								  ActiveType rho,
								  ActiveType nu )
		: S0_(S0), lambda_(lambda), alpha_(alpha), beta_(beta), rho_(rho), nu_(rho) {
			// check for valid parameter inputs
		}
		// stochastic process interface
		// dimension of X
		inline virtual size_t size() { return 2; }
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline virtual size_t factors() { return 2; }
		// initial values for simulation
		inline virtual VecP initialValues() {
			VecP X(2);
			X[0] = S0_;
			X[1] = alpha_;
			return X;
		}
		// a[t,X(t)]
		inline virtual VecA drift( const DateType t, const VecA& X) {
			VecA a(2);
			// S-variable drift-less
			a[0] = 0.0;
            // z-variable drift-less
			a[1] = 0.0;
			return a;
		}
		// b[t,X(t)]
		inline virtual MatA diffusion( const DateType t, const VecA& X) {
			MatA B(2);
			B[0].resize(2);
			B[1].resize(2);
			// S-variable z(t) (S(t) + lambda)^beta dW(t)
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

#endif  /* ifndef quantlib_templateshiftedsabrmodel_hpp */
