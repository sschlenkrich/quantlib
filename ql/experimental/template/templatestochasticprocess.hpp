/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templateprocess.hpp
    \brief define interface for general multi-dimensional stochastic process
	           
*/


#ifndef quantlib_templatestochasticprocess_hpp
#define quantlib_templatestochasticprocess_hpp


namespace QuantLib {


	// Declaration of stochastic process class
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateStochasticProcess {
	public:
		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		// subset of QL's StochasticProcess interface for X = [ x, y, z, d ] (y row-wise)
		// with dX = a[t,X(t)] dt + b[t,X(t)] dW

		// dimension of X
		inline virtual size_t size() = 0;
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline virtual size_t factors() = 0;
		// initial values for simulation
		inline virtual VecP initialValues() = 0;
		// a[t,X(t)]
		inline virtual VecA drift( const DateType t, const VecA& X) = 0;
		// b[t,X(t)]
		inline virtual MatA diffusion( const DateType t, const VecA& X) =0;

		// integrate X1 = X0 + drift()*dt + diffusion()*dW*sqrt(dt)
		// default implementation
		inline virtual void evolve( const DateType t0, const VecA& X0, const DateType dt, const VecD& dW, VecA& X1 ) {
			// ensure X1 has size of X0
			VecA a = drift(t0, X0);
			MatA b = diffusion(t0, X0);
			for (size_t i=0; i<X1.size(); ++i) {
				X1[i] = 0.0;
				for (size_t j=0; j<dW.size(); ++j) X1[i] += b[i][j]*dW[j];
				X1[i] = X0[i] + a[i]*dt + X1[i]*sqrt(dt);
			}
			return;
		}

		// default implementation, zero interest rates
		inline virtual ActiveType numeraire(const DateType t, const VecA& X) { return 1.0; }

		// default implementation, zero interest rates
		inline virtual ActiveType zeroBond(const DateType t, const DateType T, const VecA& X) { return 1.0; }

		// default implementation for single-asset models
		inline virtual ActiveType asset(const VecA& X) { return X[0]; }
		
	};

}

#endif  /* ifndef quantlib_templatestochasticprocess_hpp */
