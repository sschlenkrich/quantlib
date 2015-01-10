/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templatequasigaussian.hpp
    \brief (MC) pricing for multi-factor quasi-Gaussian model with stochastic vol
	           
			   r(t) = f(0,t) + 1^T*x(t)
			   
			   dx(t)     = [ y(t)*1 - a*x(t) ] dt                                        + sqrt[z(t)]*sigma_x^T(t,x,y) dW
			   dy(t)     = [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - a*y(t) - y(t)*a ] dt
			   dz(t)     = theta [ z0 - z(t) ] dt                                        + eta(t)*sqrt[z(t)]           dZ
			   d beta(t) = r(t)*beta(t) dt  (bank account numeraire)
			   
		   All methods are template based to allow incorporation of Automatic Differentiation
		   tools
*/


#ifndef quantlib_templatequasigaussian_hpp
#define quantlib_templatequasigaussian_hpp

#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/option.hpp>
#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>
#include <ql/experimental/template/auxilliaries/templateintegrators.hpp>
#include <ql/experimental/template/auxilliaries/templatesvd.hpp>



namespace QuantLib {

	class TemplateModel : public virtual Observable { };

	// Declaration of the quasi-Gaussian model class
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateQuasiGaussianModel : public TemplateModel {
	protected:

		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		// attributes defining the model
		Handle<YieldTermStructure> termStructure_;  // the yield curve is assumed to be passive
		// number of yield curve factors (excluding stoch. vol)
		size_t                     d_;       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
		// unique grid for time-dependent parameters
		VecD                       times_;   // time-grid of left-constant model parameter values
		// time-dependent parameters, left-piecewise constant on times_-grid
		MatA                       lambda_;  // volatility
		MatA                       alpha_;   // shift
		MatA                       b_;       // f-weighting
		VecA                       eta_;     // vol-of-vol
		// time-homogeneous parameters
		VecP                       delta_;   // maturity of benchmark rates f(t,t+delta_i) 		
		VecP                       chi_;     // mean reversions
		MatP                       Gamma_;   // (benchmark rate) correlation matrix
		// stochastic volatility process parameters
		PassiveType                theta_;   // mean reversion speed
		PassiveType                z0_;      // mean reversion level z0=z(0)=1

		// additional parameters (calculated at initialisation via SVD)
		MatP                       DfT_;     // factorized correlation matrix Df^T with Df^T * Df = Gamma
		MatP                       HHfInv_;  // weighting matrix H*Hf^-1 = [ exp{-chi_i*delta_j} ]^-1

		// lightweight container holding the current state of the yield curve
		class State {
		public:
		    VecA        x;
			MatA        y;
			ActiveType  z;
			ActiveType  beta;
			// constructor
			State( VecA X, size_t d) {
				QL_REQUIRE(X.size()==d+d*d+1+1,"TemplateQuasiGaussianModel::State Constructor: Dimensions mismatch.");
				x.resize(d);
				y.resize(d);
				for (size_t k=0; k<d; ++k) x[k] = X[k];
				for (size_t i=0; i<d; ++i) {
					y[i].resize(d);
					for (size_t j=0; j<d; ++j) y[i][j] = X[d+i*d+j];  // y row-wise
				}
				z    = X[d+d*d  ];
				beta = X[d+d*d+1];
			}
			inline void toVec(VecA& X) {
				size_t d = x.size();
				// check dimensions
				QL_REQUIRE(y.size()==d,"TemplateQuasiGaussianModel::State Assignment: y-row dimension mismatch.");
				for (size_t k=0; k<d; ++k) {
					QL_REQUIRE(y[k].size()==d,"TemplateQuasiGaussianModel::State Assignment: y-column dimension mismatch.")
				}
				X.resize(d+d*d+1+1);
				for (size_t k=0; k<d; ++k) X[k]       = x[k];
				for (size_t k=0; k<d; ++k) X[d+i*d+j] = y[i][j];
				X[d+d*d  ]                            = z;
				X[d+d*d+1]                            = beta;
			}
		};

		inline bool checkModelParameters(bool throwException=true){
			bool ok = true;
			// check yield curve...
			// non-zero dimension
			if (d_<1)               { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel number of factors larger zero required."); }
			// non-empty time-grid
			size_t n=times_.size();
			if (n<1)                { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-empty time-grid required."); }
			// dimensions of time-dependent parameters
			if (lambda_.size()!=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong lambda dimension."); }
			if (alpha_.size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong alpha dimension.");  }
			if (b_.size()     !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong b dimension.");      }
			for (size_t k=0; k<d_; ++k) {
			    if (lambda_[k].size()!=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong lambda time dimension."); }
			    if (alpha_[k].size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong alpha time dimension.");  }
			    if (b_[k].size()!=n)      { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong b time dimension.");      }
			}
			if (eta_.size() !=n) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong eta time dimension.");    }
			// dimensions of time-homogeneous parameters
			if (delta_.size()!=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong delta dimension."); }
			if (chi_.size()  !=d_)   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong chi i-dimension.");   }
			if (Gamma_.size()  !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma j-dimension."); }
			for (size_t k=0; k<d_; ++k) {
			    if (Gamma_[k].size() !=d_) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel wrong Gamma dimension."); }
			}
			// plausible parameter values
			// ascending time-grid
			for (size_t k=0; k<times_.size()-1; ++k) {
				if (times_[k]>=times_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending time-grid required."); }
			}
			// non-negative values
			for (size_t j=0; j<n; ++j) {
				for (size_t i=0; i<d_; ++i) {
                    if (lambda_[i][j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel lambda>=0 required."); }
                    if (alpha_[i][j] <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel alpha>=0 required.");  }
                    if (b_[i][j]     <0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel b>=0 required.");      }
				}
                if (eta_[j]<0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel eta>=0 required.");      }
			}
			// positive/ascending values
			if (delta_[0]<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel delta>0 required."); }
			if (chi_[0]  <=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel chi>0 required."); }
			for (size_t k=0; k<d_-1; ++k) {
				if (delta_[k]>=delta_[k+1]) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending delta values required."); }
				if (chi_[k]  >=chi_[k+1])   { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel ascending chi values required."); }
			}
			// plausible correlation values
			for (size_t i=0; i<d_; ++i) {
				for (size_t j=i; j<d_; ++j) {
					if (Gamma_[i][j]<0.0) {
						ok = false;
						if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][j]>=0 required."); 
					}
					if (i==j) {
						if (Gamma_[i][j]!=1.0) {
							ok = false;
							if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][i]=1 required."); 
						}
					}
					if (Gamma_[i][j]!=Gamma_[j][i]) {
						ok = false;
						if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma[i][j]=Gamma[j][i] required."); 
					}
					if (i<j) {
						if((Gamma_[i][j-1]>=Gamma_[i][j])|(Gamma_[i+1][j]>=Gamma_[i][j])) {
							ok = false;
							if (throwException) QL_REQUIRE(false,"QuasiGaussianModel Gamma descending sub-diagonals required."); 
						}
					}
				}
			}
			// stochastic vol parameters
			if (theta_<=0.0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel theta>0 required."); }
			if (z0_!=1.0)    { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel z0=1 required.");    }
            // finished
			return ok;
		}

		// evaluate Df^T with Df^T * Df = Gamma and H*Hf^-1 via singular value decomposition
		// return false (and throw exception) on error
		inline bool factorMatrices(bool throwException=true){
			bool ok = true;
			// row-wise matrices
			size_t dim = d_;
			PassiveType *A  = new PassiveType[dim*dim];
			PassiveType *U  = new PassiveType[dim*dim];
			PassiveType *S  = new PassiveType[dim*dim];
			PassiveType *VT = new PassiveType[dim*dim];
			// dummy auxilliary variables
			PassiveType work;
			int lwork, info;
			// Gamma = U S V^T
			for (size_t i=0; i<dim; ++i) {
				for (size_t j=0; j<dim; ++j) {
					A[i*dim+j] = Gamma_[i][j];
				}
			}
			TemplateAuxilliaries::svd("S","S",(int*)&dim,(int*)&dim,A,(int*)&dim,S,U,(int*)&dim,VT,(int*)&dim,&work,&lwork,&info);
			// check min(S)>0
			PassiveType minS=S[0];
			for (size_t i=1; i<dim; ++i) if (S[i*dim+i]<minS) minS = S[i*dim+i];
			if (minS<=0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-singular Gamma required."); }
			// evaluate Df^T = U S^{1/2}
			DfT_.resize(dim);
			for (size_t i=0; i<dim; ++i) {
				DfT_[i].resize(dim);
				for (size_t j=0; j<dim; ++j) {
					DfT_[i][j] = U[i*dim+j]*sqrt(S[j*dim+j]);
				}
			}
			// H*Hf = [ exp{-chi_i*delta_j} ] = U S V^T
			for (size_t i=0; i<dim; ++i) {
				for (size_t j=0; j<dim; ++j) {
					A[i*dim+j] = exp(-chi_[i]*delta_[j]);
				}
			}
			TemplateAuxilliaries::svd("S","S",(int*)&dim,(int*)&dim,A,(int*)&dim,S,U,(int*)&dim,VT,(int*)&dim,&work,&lwork,&info);
			// check min(S)>0
			minS=S[0];
			for (size_t i=1; i<dim; ++i) if (S[i*dim+i]<minS) minS = S[i*dim+i];
			if (minS<=0) { ok = false; if (throwException) QL_REQUIRE(false,"QuasiGaussianModel non-singular Gamma required."); }
			// evaluate H*Hf^-1 = V S^{-1} U^T
			for (size_t i=0; i<dim; ++i) {
				for (size_t j=0; j<dim; ++j) {
					HHfInv_[i][j] = 0.0;
					for (size_t k=0; k<dim; ++k) HHfInv_[i][j] += U[j*dim+k] * VT[k*dim+j] / S[k*dim+k];
				}
			}
			// finished
			return ok;
		}

		public:  // for debugging purpose we allow unsafe aaccess to restricted members, IN GENERAL NO CHECK FOR DIMENSIONS!

		// Constructor
		TemplateQuasiGaussianModel() { checkModelParameters(); factorMatrices(); }

		TemplateQuasiGaussianModel(
			const Handle<YieldTermStructure>& termStructure,
		    // number of yield curve factors (excluding stoch. vol)
		    size_t                     d,       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
		    // unique grid for time-dependent parameters
		    const VecD &                times_,   // time-grid of left-constant model parameter values
		    // time-dependent parameters, left-piecewise constant on times_-grid
		    const MatA &                lambda,  // volatility
		    const MatA &                alpha,   // shift
		    const MatA &                b,       // f-weighting
		    const VecA &                eta,     // vol-of-vol
		    // time-homogeneous parameters
		    const VecP &                delta,   // maturity of benchmark rates f(t,t+delta_i) 		
		    const VecP &                chi,     // mean reversions
		    const MatP &                Gamma,   // (benchmark rate) correlation matrix
		    // stochastic volatility process parameters
		    PassiveType                theta   // mean reversion speed
			) : termStructure_(termStructure), d_(d), times_(times), lambda_(lambda), alpha_(alpha), b_(b), eta_(eta),
			    delta_(delta), chi_(chi), Gamma_(Gamma), theta_(theta), z0_((PassiveType)1.0) {
                checkModelParameters();
				// calculate  DfT_
				// calculate  HHfInv_
				factorMatrices();
			}

		// helpers

		inline size_t maxidx( size_t i ) { return (i<d_) ? i : d_-1; }

	    // evaluate n s.t. t[n-1] < t <= t[n]
		inline size_t idx( DateType t) {
			if ((t <= times_[0]) | (times_.size()<2)) return 0;
			if (t >  times_[times_.size()-2 ])        return times_.size()-1;
			// bisection search
			size_t a = 0, b = times_.size()-2
			while (b-a>1) {
			    size_t s = (a + b) / 2;
				if (t <= times_[s]) b = s;
				else                a = s;
			}
			return b;
		}

		// parameter functions (no dimension checks)
		inline ActiveType lambda( size_t i, DateType t) { return lambda_[maxidx(i)][idx(t)]; }
		inline ActiveType alpha ( size_t i, DateType t) { return alpha_[maxidx(i)][idx(t)];  }
		inline ActiveType b     ( size_t i, DateType t) { return b_[maxidx(i)][idx(t)];      }
		inline ActiveType eta   ( DateType t)           { return eta_[idx(t)];       }

		// analytic formulas

		inline ActiveType G(size_t i, DateType t, DateType T) { return (1.0-exp(-chi_[i]*(T-t)))/chi_[i]; }

		inline
        ActiveType shortRate ( DateType t, const VecA& x ) {
		    ActiveType r = termStructure_forwardRate(t,t,Continuous);
			for (size_t k=0; k<d_; ++k) r += x[k];
			return r;
        }

		inline
		ActiveType forwardRate( DateType t, DateType T, const VecA& x, const MatA&  y) {
			ActiveType f = termStructure_forwardRate(t,T,Continuous);
			for (size_t i=0; i<d_; ++i) {
				ActiveType tmp = x[i];
				for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
				f += exp(-chi_[i]*(T-t)) * tmp;
			}
			return f;
		}

		inline
		ActiveType ZeroBond( DateType t, DateType T, const VecA& x, const MatA&  y) {
		    PassiveType DF1  = termStructure_->discount(t);
		    PassiveType DF2  = termStructure_->discount(T);
		    ActiveType  Gx   = 0;   // G^T * x
			for (size_t i=0; i<d_; ++i) Gx += x[i]*G(i,t,T);
            ActiveType  GyG  = 0;   // G^T * y * G
			for (size_t i=0; i<d_; ++i) {
				ActiveType tmp = 0;
				for (size_t j=0; j<d_; ++j) tmp += y[i][j]*G(j,t,T);
				GyG += G(i,t,T)*tmp;
			}
			ActiveType ZCB = DF2 / DF1 * exp(-Gx - 0.5*GyG);
			return ZCB;
		}

		inline  // diagonal vector
		VecA sigma_f( DateType t, const VecA& x, const MatA&  y) {
			VecA res(d_);
			for (size_t k=0; k<d_; ++k) res[k] = lambda(k,t) * (alpha(k,t) + b(k,t)*forwardRate(t,t+delta_[k],x,y));
			return res;
		}

		inline  // sigma_x^T
		MatA sigma_xT( DateType t, const VecA& x, const MatA&  y) {
			MatA tmp(d_), res(d_);
			VecA sigmaf = sigma_f(t,x,y);
			// tmp = sigma_f * Df^T
			for (size_z i=0; i<d_; ++i) {
				tmp[i].resize(d_);
				for (size_t j=0; j<d_; ++j) {
					tmp[i][j] = sigmaf[i] * DfT_[i][j];
				}
			}
			// res = H*Hf^-1 * tmp
			for (size_z i=0; i<d_; ++i) {
				res[i].resize(d_);
				for (size_t j=0; j<d_; ++j) {
					res[i][j] = 0;
					for (size_t k=0; k<d_; ++k) res[i][j] += HHfInv_[i][k]*tmp[k][j];
				}
			}
			return res;
		}


		// subset of QL's StochasticProcess interface for X = [ x, y, z, d ] (y row-wise)
		// with dX = a[t,X(t)] dt + b[t,X(t)] dW

		// dimension of X
		inline size_t size()    { return d_ + d_*d_ + 1 + 1; }
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline size_t factors() { return d_ + 1; }
		// initial values for simulation
		inline VecP initialValues() {
			VecP X(size());
			for (size_t k=0; k<d_ + d_*d_; ++k)       X[k] = 0.0;  // x(0), y(0)
			for (size_t k=d_ + d_*d_; k<size(); ++k)  X[k] = 1.0;  // z(0), beta(0)
			return X;
		}

		// a[t,X(t)]
		inline VecA drift( DateType t, VecA X) {
			VecA a(size());
			State state(X,d_);
			// x-variable [ y(t)*1 - chi*x(t) ]
			for (size_t k=0; k<d_; ++k) {
				a[k] = -chi_[k]*state.x[k];
				for (size_t j=0; j<d_; ++j) a[k] += state.y[k][j];
			}
			// y-variable [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - chi*y(t) - y(t)*chi ]
			MatA sigmaxT = sigma_xT(t,state.x,state.y);
			for (i=0; i<d_; ++i) {
				for (j=0; j<d_; ++j) {
					a[d_+i*d_+j] = 0.0;
					for (size_t k=0; k<d_; ++k) a[d_+i*d_+j] += sigmaxT[i][k]*sigmaxT[k][j];
					a[d_+i*d_+j] *= state.z;
					a[d_+i*d_+j] -= (chi_[i]+chi_[j])*state.y[i][j];
				}
			}
			// z-variable theta [ z0 - z(t) ]
			a[d_+d_*d_] = theta_*(z0_ - state.z);
			// beta-variable r(t)*beta(t)
			a[d_+d_*d_+1] = shortRate(t,state.x)*state.beta;
			// finished
			return a;
		}

		// b[t,X(t)]
		inline MatA diffusion( DateType t, VecA X) {
			Mat b(size());
			for (size_t k=0; k<size(); ++k) b[k].resize(factors());
			State state(X,d_);
			ActiveType sqrtz = sqrt(abs(state.z));
			MatA sigmaxT = sigma_xT(t,state.x,state.y);
			// x-variable sqrt[z(t)]*sigma_x^T(t,x,y)
			for (size_t i=0; i<d_; ++i) {
				for (size_t j=0; j<d_; ++j) b[i][j] =  sqrtz * sigmaxT[i][j];
				b[i][d] = 0.0;
			}
			// y-variable 0
			for (size_t i=d_; i<d_+d_*d_; ++i) {
				for (size_t j=0; j<d_+1; ++j) {
					b[i][j] = 0.0;
				}
			}
			// z-variable eta(t)*sqrt[z(t)]
			for (size_t j=0; j<d_+1; ++j) b[d_+d_*d_][j] = 0.0;
			b[d_+d_*d_][d] = eta(t)*sqrtz;
			// beta-variable 0
			for (size_t j=0; j<d_+2; ++j) b[d_+d_*d_][j] = 0.0;
			// finished
			return b;
		}


	};

}

#endif  /* ifndef quantlib_templatequasigaussian_hpp */
