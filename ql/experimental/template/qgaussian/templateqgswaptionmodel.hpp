/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file templateqgswaptionmodel.hpp
    \brief (Approximate) swaption pricing for multi-factor quasi-Gaussian model with stochastic vol
	           
			   r(t) = f(0,t) + 1^T*x(t)
			   
			   dx(t)     = [ y(t)*1 - a*x(t) ] dt                                        + sqrt[z(t)]*sigma_x^T(t,x,y) dW
			   dy(t)     = [ z(t)*sigma_x^T(t,x,y)*sigma_x(t,x,y) - a*y(t) - y(t)*a ] dt
			   dz(t)     = theta [ z0 - z(t) ] dt                                        + eta(t)*sqrt[z(t)]           dZ
			   ds(t)     = r(t) dt  ( s(t) = int_0^t r(s) ds, for bank account numeraire)
			   
		   All methods are template based to allow incorporation of Automatic Differentiation
		   tools
*/


#ifndef quantlib_templateqgswaptionmodel_hpp
#define quantlib_templateqgswaptionmodel_hpp

#include <boost/shared_ptr.hpp>
#include <ql/types.hpp>

#include <ql/experimental/template/qgaussian/templatequasigaussian.hpp>
//#include <ql/experimental/template/auxilliaries/gausslobatto.hpp>



namespace QuantLib {

	// Declaration of the quasi-Gaussian model class
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateQGSwaptionModel {
	protected:

		// container class definitions
		typedef std::vector<DateType>                      VecD;
		typedef std::vector<PassiveType>                   VecP; 
		typedef std::vector<ActiveType>                    VecA;
		typedef std::vector< std::vector<DateType> >       MatD;
		typedef std::vector< std::vector<PassiveType> >    MatP;
		typedef std::vector< std::vector<ActiveType> >     MatA;

		// reference to QG modelspecs
		boost::shared_ptr< TemplateQuasiGaussianModel<DateType,PassiveType,ActiveType> >   model_;

		// time grid for barX, barY evaluation
		VecD               times_;   
		MatA               barX_;     // E^A [ x(t) ] 
		std::vector<MatA>  barY_;     // E^A [ y(t) ] 
		VecA               x0_;       // = VecA(model_->factors()-1,0.0);
		MatA               y0_;       // = MatA(model_->factors()-1,x0_);


		// swap spec container
		struct Swap {
			VecD   times;    // T[0], ..., T[N]
			VecD   weights;  // w[0], ..., w[N-1]
			size_t N;
			Swap ( const VecD& t, const VecD& w) : times(t), weights(w), N(weights.size()) {}
		};

		// E^A [ y(T) ] = H(T) H(t)^-1 E^A [ y(t) ] H(t)^-1 H(T) + int_t^T H(T) H(s)^-1 sigma_x^T sigma_x H(s)^-1 H(T) ds
		//              = H(T) H(t)^-1 E^A [ y(t) ] H(t)^-1 H(T) + int_t^T [sigma_x^T sigma_x]_i,j exp{-(chi_i + chi_j)(t-s)} |_i,j=1,..,d
		// via mid-point rule
		inline MatA expectationAy( DateType t, DateType T, const MatA& barYt ) {
			MatA barYT = barYt;
			// H(T) H(t)^-1 E^A [ y(T) ] H(t)^-1 H(T)
			for (size_t i=0; i<model_->factors()-1; ++i) {
       			for (size_t j=0; j<model_->factors()-1; ++j) {
					barYT *= exp(-(model_->chi()[i] + model_->chi()[j]) * (T-t) );
				}
			}
			// int_t^T H(T) H(s)^-1 sigma_x^T sigma_x H(s)^-1 H(T) ds
			MatA sig_xT = model_->sigma_xT((t+T)/2.0,x0_,y0_);
			for (size_t i=0; i<model_->factors()-1; ++i) {
       			for (size_t j=0; j<model_->factors()-1; ++j) {
					// integrand at (t+T)/s
                    ActiveType tmp = 0.0;
					for (size_t k=0; k<model_->factors()-1; ++k) tmp += sig_xT[i][k] * sig_xT[j][k];
                    tmp *= exp( -(model_->chi()[i] + model_->chi()[j]) * (T-t)/2.0 );
					// mid-point rule
					barYT[i][j] += tmp * (T - t);
				}
			}
			return MatYT;
		}

		// E^A [ x(T) ] = H(T)H(t)^-1 E^A [ x(t) ] + int_t^T H(T)H(s)^-1 [ E^A[y(s)]*1 + sigma_x^T sigma_x G_A(s) ] ds
		inline VecA expectationAx( const Swap& s, DateType t, DateType T, const MatA& barYtT, const VecA& barXt ) {
			VecA barXT = barXt;
			// H(T)H(t)^-1 E^A [ x(t) ]
			for (size_t i=0; i<model_->factors()-1; ++i) barXT[i] *= exp(-model_->chi()[i] * (T-t) );
			// int_t^T H(T)H(s)^-1 [ E^A[y(s)]*1 + sigma_x^T sigma_x G_A(s) ] ds
			// G_A((t+T)/2)
			VecA GA(model_->factors()-1);
			for (size_t i=0; i<model_->factors()-1; ++i) {
				GA[i] = 0.0;
				for (size_t j=0; j<s.N; ++j) GA[i] += s.weights[j] * model_->zeroBond(0.0,s.times[j+1],x0_,y0_) * model_->G(i,(t+T)/2.0,s.times[j+1]);
				GA[i] *= 1.0 / annuity(s,0.0,x0_,y0_);
			}
			// v = sigma_x^T [sigma_x G_A(s)]  | s=(t+T)/2
			MatA sig_xT = model_->sigma_xT((t+T)/2.0,x0_,y0_);
			VecA u(model_->factors()-1);
			for (size_t i=0; i<model_->factors()-1; ++i) {
				u[i] = 0.0;
				for (size_t j=0; j<model_->factors()-1; ++j) u[i] += sig_xT[j][i]*GA[j];
			}
			VecA v(model_->factors()-1);
			for (size_t i=0; i<model_->factors()-1; ++i) {
				v[i] = 0.0;
				for (size_t j=0; j<model_->factors()-1; ++j) v[i] += sig_xT[i][j]*u[j];
				// E^A[y(s)]*1 + v  | s=(t+T)/2
				for (size_t j=0; j<model_->factors()-1; ++j) v[i] += barYtT[i][j];
				// H(T)H(s)^-1 E^A[y(s)]*1 + v  | s=(t+T)/2
				v[i] *= exp( -model_->chi()[i] * (T-t)/2.0 );
				// mid-point rule
				barXT[i] += v[i] * (T - t);
			}
			return barXT;
		}

		// annuity
		inline ActiveType annuity( const Swap& s, DateType t, const VecA& x, const MatA& y) {
			ActiveType den = 0.0;
			for (size_t k=0; k<s.N; ++k) den += s.weights[k]*model_->zeroBond(t,s.times[k+1],x,y);
		}

		// forward swaprate
		inline ActiveType swapRate( const Swap& s, DateType t, const VecA& x, const MatA& y) {
			ActiveType num = model_->zeroBond(t,s.times[0],x,y) - model_->zeroBond(t,s.times[s.N],x,y);
			return num / annuity(s,t,x,y);
		}

		// d P(t,T,x,y) / dx  = -P(t,T,x,y) G(t,T) 
		inline VecA zcbGradient( DateType t, DateType T, const VecA& x, const MatA& y) {
			VecA grad(x.size());
			ActiveType DF = model_->zeroBond(t,T,x,y);
			for (size_t k=0; k<grad.size(); ++k) grad[k] = -DF * model_->G(k,t,T);
			return grad;
		}

		// gradient 
		VecA swapGradient( const Swap& s, DateType t, const VecA& x, const MatA& y) {
			VecA grad(x.size());
			VecA gZCB(x.size());
			// numerator
			grad = zcbGradient(t,s.times[0],x,y);
			gZCB = zcbGradient(t,s.times[s.N],x,y);
			ActiveType annuit = annuity(s,t,x,y);
			for (size_t k=0; k<grad.size(); ++k) grad[k] = (grad[k] - gZCB[k]) / annuit;
			// denumerator
			ActiveType dSdAnnuity = - swapRate(s,t,x,y) / annuit;
			for (size_t i=0; i<s.N; ++i) {
				gZCB = zcbGradient(t,s.times[i+1],x,y);
				for (size_t k=0; k<grad.size(); ++k) grad[k] += dSdAnnuity * s.weights[i] * gZCB[k];
			}
			return grad;
		}

		// hessian via central differences to avoid errors in tedious calculations
		MatA swapHessian( const Swap& s, DateType t, const VecA& x, const MatA& y) {
			PassiveType eps = 1.0-7;
			MatA hess(x.size());
			VecA xm(x), xp(x), ym(x.size()), yp(x.size());
			for (size_t k=0; k<x.size(); ++k) {
				// bump
				xm[k] -= eps;
				xp[k] += eps;
				// recalculate
				ym = swapGradient(s,t,xm,y);
				yp = swapGradient(s,t,xp,y);
				// second derivs
				hess[k].resize(x.size());
				for (size_t i=0; i<x.size(); ++i) hess[k][i] = yp[i]/2.0/eps - ym[i]/2.0/eps;
				// reset
				xm[k] = x[k];
				xp[k] = x[k];
			}
			// maybe better check symetry
			return hess;
		}

		// lambda_S for given x,y
		ActiveType lambda( const Swap& s, DateType t, const VecA& x, const MatA& y) {
			VecA grad  = swapGradient(s,t,x,y);
			MatA sigxT = model_->sigma_xT(t,x,y);
			VecA u(x.size());
			// grad * sigma_x^T
			for (size_t j=0; j<x.size(); ++j) {
				u[j] = 0.0;
				for (size_t i=0; i<x.size(); ++i) u[j] += grad[i] * sigxT[i][j];
			}
			// [ grad sigma_x^T sigma_x grad^T ]^{1/2}
			ActiveType num = 0.0;
			for (size_t j=0; j<x.size(); ++j) num += u[j]*u[j];
			num = sqrt(num);
			// swapt rate at (t,x,y)=0
			VecA x0(x.size(),(ActiveType)0.0);
			MatA y0(x.size());
			for (size_t i=0; i<x.size(); ++i) y0[i].resize(x.size(),(ActiveType)0.0);
			ActiveType den = swapRate(s,0.0,x0,y0);
			return num / den;
		}

		// b_S for given x,y
		ActiveType b( const Swap& s, DateType t, const VecA& x, const MatA& y) {

			// swapt rate at (t,x,y)=0
			VecA x0(x.size(),(ActiveType)0.0);
			MatA y0(x.size());
			for (size_t i=0; i<x.size(); ++i) y0[i].resize(x.size(),(ActiveType)0.0);
			ActiveType S0 = swapRate(s,0.0,x0,y0);

			// u = grad * sigma_x^T
			VecA u(x.size());
			VecA grad  = swapGradient(s,t,x,y);
			MatA sigxT = model_->sigma_xT(t,x,y);
			for (size_t j=0; j<x.size(); ++j) {
				u[j] = 0.0;
				for (size_t i=0; i<x.size(); ++i) u[j] += grad[i] * sigxT[i][j];
			}
			// den = [ grad sigma_x^T sigma_x grad^T ]^2
			ActiveType den = 0.0;
			for (size_t i=0; i<x.size(); ++i) den += u[i]*u[i];
			den = den * den;

			// d_x ...

			// v = grad c_x = grad sigma_x^T sigma_x = u * sigma_x = sigma_x^T u^T
			VecA v(x.size());
			for (size_t i=0; i<x.size(); ++i) {
				v[i] = 0.0;
				for (size_t j=0; j<x.size(); ++j) v[i] += sigxT[i][j] * u[j];
			}

			// inputs...
			MatP HHfInv = model->HHfInv();
			MatP DfT    = model->DfT();
			VecP delta  = model->delta();
			VecP chi    = model->chi();
			VecA lambda_f(x.size());
			VecA b_f(x.size());
			for (size_t k=0; k<x.size(); ++k) { 
				lambda_f[k] = model->lambda(k,t);
				b_f[k]      = model->b(k,t);
			}

			// d sigma_x^T / d x_k
			MatA dx(x.size());
			MatA dsigmaxT_dxk(x.size());
			VecA diag(x.size());
			MatA tmp(x.size());
			for (size_t k=0; k<x.size(); ++k) {
				dx[k].resize(x.size(),(ActiveType)0.0);
				dsigmaxT_dxk[k].resize(x.size());
				tmp[k].resize(x.size());
			}

            // d_x = sum...
			for (size_t k=0; k<x.size(); ++k) {
			    // diag^k = diag_i [ lambda_{f,i} b_{f,i} h_k(t+delta_i)/h_k(t) ]
			    //        = diag_i [ lambda_{f,i} b_{f,i} exp(-chi_k delta_i)   ]
			    for (size_t i=0; i<x.size(); ++i) diag[i] = lambda_f[i] * b_f[i] * exp(-chi[k]*delta[i]);
				// tmp = diag * DfT
				for (size_t i=0; i<x.size(); ++i) {
					for (size_t j=0; j<x.size(); ++j) tmp[i][j] = diag[i] * DfT[i][j];
				}
				// d sigma_x^T / d x_k = HHf^-1 * diag * DfT = HHf^-1 * tmp
				for (size_t i=0; i<x.size(); ++i) {
					for (size_t j=0; j<x.size(); ++j) {
						dsigmaxT_dxk[i][j] = 0.0;
						for (size_t l=0; l<x.size(); ++l) dsigmaxT_dxk[i][j] += HHfInv[i][l] * tmp[l][j];
					}
				}
				// tmp = [d sigma_x^T / d x_k] sigma_x = [d sigma_x^T / d x_k] [sigma_x^T]^T
				for (size_t i=0; i<x.size(); ++i) {
					for (size_t j=0; j<x.size(); ++j) {
						tmp[i][j] = 0.0;
						for (size_t l=0; l<x.size(); ++l) tmp[i][j] += dsigmaxT_dxk[i][l] * sigxT[j][l];
					}
				}
				// d_x = sum_k v_k [tmp + tmp^T]
				for (size_t i=0; i<x.size(); ++i) {
					for (size_t j=0; j<x.size(); ++j) dx[i][j] += v[k] * (tmp[i][j] + tmp[j][i]); 
				}

			}  // for k=0..d

			// num1 and num2...

			// num1 = grad * (dx grad^T)  (overwriting vector u)
			for (size_t i=0; i<x.size(); ++i) {
				u[i] = 0.0;
				for (size_t j=0; j<x.size(); ++j) u[i] += dx[i][j] * grad[j];
			}
			ActiveType num1 = 0.0
			for (size_t i=0; i<x.size(); ++i)  num1 += grad[i] * u[i];

			// num2 = grad c_x grad^2 c_x grad^T = v (grad^2 v^T)
			tmp = swapHessian(s,t,x,y);
			for (size_t i=0; i<x.size(); ++i) {
				u[i] = 0.0;
				for (size_t j=0; j<x.size(); ++j) u[i] += tmp[i][j] * v[j];
			}
			ActiveType num2 = 0.0
			for (size_t i=0; i<x.size(); ++i)  num2 += v[i] * u[i];

			// finishing up
			ActiveType res = S0 * ( num1/den/2.0 + num2/den );
			return res;
		}

	public:
		TemplateQGSwaptionModel (
			const boost::shared_ptr< TemplateQuasiGaussianModel<DateType,PassiveType,ActiveType> >&   model )
			: model_(model) { }


	};


}

#endif  /* ifndef quantlib_templateqgswaptionmodel_hpp */
