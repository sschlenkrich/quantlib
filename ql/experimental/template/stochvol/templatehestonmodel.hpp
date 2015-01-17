/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_templatehestonmodels_hpp
#define quantlib_templatehestonmodels_hpp

#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <ql/errors.hpp>
#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>
#include <ql/experimental/template/auxilliaries/gausslobatto.hpp>
#include <ql/experimental/template/auxilliaries/Complex.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    //using namespace std;
	
	// heston model with analytic vanilla pricing 
	//    dS(t) = S(t) sqrt[ v(t) ] dW(t)
    //    dv(t) = kappa [theta - v(t)] dt + sigma sqrt[ v(t) ] dZ(t)
    //    dW(t) dZ(t) = rho dt

	template <class DateType, class PassiveType, class ActiveType>
	class TemplateHestonModel {
    //typedef std::complex<Real> complex;
    typedef Cpx::Complex<ActiveType> complex;
	protected:
        ActiveType kappa_;    // mean reversion speed of volatility
        ActiveType theta_;    // mean reversion level of volatility
        ActiveType sigma_;    // volatility of volatility
        ActiveType rho_;      // correlation vol vs. underlying
        ActiveType v0_;       // initial volatility
    public:
        // constructor
        TemplateHestonModel( ActiveType kappa,
                             ActiveType theta,
                             ActiveType sigma,
                             ActiveType rho,
                             ActiveType v0 ) :
        kappa_(kappa), theta_(theta), sigma_(sigma), rho_(rho), v0_(v0) {}
                              
        // inspectors
        const ActiveType& kappa() const { return kappa_; }
        const ActiveType& theta() const { return theta_; }
        const ActiveType& sigma() const { return sigma_; }
        const ActiveType& rho()   const { return rho_;   }
        const ActiveType& v0()    const { return v0_;    }
        // maths
        bool fellerConstraint() {
            return (sigma >= 0.0 && sigma*sigma < 2.0*kappa*theta);
        }
        // undiscounted expectation of vanilla payoff
        ActiveType vanillaOption(const PassiveType forwardPrice,
                                 const PassiveType strikePrice,
                                 const DateType    term,
                                 const int         callOrPut,
                                 const PassiveType accuracy,
                                 const size_t      maxEvaluations) {
            const ActiveType c_inf = _MIN_(10.0, _MAX_(0.0001, sqrt(1.0-rho_*rho_)/sigma_)) * (v0_ + kappa_*theta_*term);
            TemplateAuxilliaries::GaussLobatto<ActiveType> integrator(maxEvaluations, accuracy);
            IntegrandGatheral gatheral1( *this, forwardPrice, strikePrice, term, 1 );
            IntegrandGatheral gatheral2( *this, forwardPrice, strikePrice, term, 2 );
            IntegrandTransformation transf1( c_inf, gatheral1 );
            IntegrandTransformation transf2( c_inf, gatheral2 );
            const ActiveType p1 = integrator.integrate( transf1, 0, 1) / M_PI; 
            const ActiveType p2 = integrator.integrate( transf2, 0, 1) / M_PI; 
            switch (callOrPut) {
                case +1:  // Call
                    return forwardPrice*(p1+0.5) - strikePrice*(p2+0.5);
                    break;
                case -1:  // Put
                    return forwardPrice*(p1-0.5) - strikePrice*(p2-0.5);
                    break;
                default:
                    QL_FAIL("unknown option type");
            }
        }

        // f_1/2( u(x) ) / ( x c_inf ), u(x) = -ln(x) / c_inf
        class IntegrandTransformation  {
        protected:
            ActiveType                                c_inf_;
            boost::function<ActiveType (ActiveType)>  f_12_;
        public:
            IntegrandTransformation(  const ActiveType& c_inf, const boost::function<ActiveType (ActiveType)>& f_12)
                : f_12_(f_12), c_inf_(c_inf) { }
            ActiveType operator()(ActiveType x) const {
                if (x * c_inf_ < QL_EPSILON) return 0;
                else                         return f_12_( -log(x) / c_inf_ ) / x / c_inf_;
            }

        };

        // key issue of Heston model, implements f_1/2 (phi) according to Gatherals approach
        class IntegrandGatheral {
        protected:
            Size j_;                                 // evaluate f_1 or f_2
            const ActiveType kappa_, theta_, sigma_, v0_;  // copy of model parameters
            // helper variables
            const DateType   term_;
            const ActiveType x_, sx_, dd_;
            const ActiveType sigma2_, rsigma_;
            const ActiveType t0_;
        public:
            // constructor
            IntegrandGatheral( 
                const TemplateHestonModel<DateType,PassiveType,ActiveType> &model,
                const ActiveType forward,     // underlying initial state
                const ActiveType strike,      // vanilla option strike
                const DateType   term,        // time to maturity
                const Size       j            // f1 or f2
                ) :
                kappa_(model.kappa()), theta_(model.theta()), sigma_(model.sigma()), v0_(model.v0()),
                term_(term),
                j_(j),
                x_(log(forward)),
                sx_(log(strike)),
                dd_(x_),
                sigma2_(sigma_*sigma_),
                rsigma_(model.rho()*sigma_),
                t0_(kappa_ - ((j== 1)? rsigma_ : 0)) {  }
            // ...
            ActiveType operator()(ActiveType phi) const {
                const ActiveType rpsig(rsigma_*phi);
                const complex t1 = t0_+complex(0, -rpsig);
                const complex d  = sqrt(t1*t1 - sigma2_*phi*complex(-phi, (j_== 1)? 1 : -1));
                const complex ex = exp(-d*ActiveType(term_));
                const complex addOnTerm =  0.0;
                if (phi != 0.0) {
                    if (sigma_ > 1e-5) {
                        const complex p = (t1-d)/(t1+d);
                        const complex g = log((1.0 - p*ex)/(1.0 - p));
                        return exp(v0_*(t1-d)*(1.0-ex)/(sigma2_*(1.0-ex*p))
                                    + (kappa_*theta_)/sigma2_*((t1-d)*term_-2.0*g)
                                    + complex(0.0, phi*(dd_-sx_))
                                    + addOnTerm
                                  ).imag()/phi;
                    }
                    else {
                        const complex td = phi/(2.0*t1)*complex(-phi, (j_== 1)? 1 : -1);
                        const complex p  = td*sigma2_/(t1+d);
                        const complex g  = p*(1.0-ex);
                        return exp(v0_*td*(1.0-ex)/(1.0-p*ex)
                                         + (kappa_*theta_)*(td*term_-2.0*g/sigma2_)
                                         + complex(0.0, phi*(dd_-sx_))
                                         + addOnTerm
                                       ).imag()/phi;
                    }
                }
                else {
                    // use l'Hospital's rule to get lim_{phi->0}
                    if (j_ == 1) {
                        const ActiveType kmr = rsigma_-kappa_;
                        if (fabs(kmr) > 1e-7) {
                            return dd_-sx_ + (exp(kmr*term_)*kappa_*theta_
                                   -kappa_*theta_*(kmr*term_+1.0) ) / (2*kmr*kmr)
                                   - v0_*(1.0-exp(kmr*term_)) / (2.0*kmr);
                        }
                        else
                            // \kappa = \rho * \sigma
                            return dd_-sx_ + 0.25*kappa_*theta_*term_*term_ + 0.5*v0_*term_;
                    }
                    else {
                        return dd_-sx_ - (exp(-kappa_*term_)*kappa_*theta_
                               + kappa_*theta_*(kappa_*term_-1.0))/(2*kappa_*kappa_)
                               - v0_*(1.0-exp(-kappa_*term_))/(2*kappa_);
                    }
                }        
                return 0;
            }  // operator()

        }; // class IntegrandGatheral

	};  // TemplateHestonModel        

	// general stochastic volatility model with constant parameters and analytic vanilla pricing formula
	//
	//    dS(t) = lambda [ b S(t) + (1-b) L ] sqrt[z(t)] dW(t)
	//    dz(t) = theta [ m - z(t) ] dt + eta sqrt[z(t)] dZ(t)
	//    dW(t) dZ(t) = rho dt
	//
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateStochVolModel {
	protected:
		boost::shared_ptr< TemplateHestonModel<DateType,PassiveType,ActiveType> > hestonModel_;
		ActiveType                                                                shift_;
	public:
		TemplateStochVolModel ( const ActiveType   lambda,
			                    const ActiveType   b,
								const ActiveType   L,
								const ActiveType   theta,
								const ActiveType   m,
								const ActiveType   eta,
								const ActiveType   z0,
								const ActiveType   rho )
		: hestonModel_(new TemplateHestonModel<DateType,PassiveType,ActiveType>(
		    // state transformations ~S(t) = S(t) + (1-b)/b L, v(t) = z(t) lambda^2 b^2
			theta,                // kappa
			m*lambda*lambda*b*b,  // theta
			eta*lambda*b,         // sigma
			rho,                  // rho
			z0*lambda*lambda*b*b  // v0
			) ), shift_( (1.0-b)/b*L ) {}

        // undiscounted expectation of vanilla payoff
        ActiveType vanillaOption(const PassiveType forwardPrice,
                                 const PassiveType strikePrice,
                                 const DateType    term,
                                 const int         callOrPut,
                                 const PassiveType accuracy,
                                 const size_t      maxEvaluations) {
			return hestonModel_->vanillaOption( forwardPrice+shift_, strikePrice+shift_, term, callOrPut, accuracy, maxEvaluations );
		}
	};


    // general stochastic volatility model interface with time-dependent parameters and parameter averaging pricing formula
	//
	//    dS(t) = lambda(t) [ b(t) S(t) + (1-b(t)) L ] sqrt[z(t)] dW(t)
	//    dz(t) = theta [ m - z(t) ] dt + eta(t) sqrt[z(t)] dZ(t)
	//    dW(t) dZ(t) = rho dt
	//
	template <class DateType, class PassiveType, class ActiveType>
	class TemplateTimeDependentStochVolModel {
	public:
		// inspectors
        virtual ActiveType  lambda( const DateType t) = 0;
	    virtual ActiveType  b(      const DateType t) = 0;
	    virtual ActiveType  L()                       = 0;
	    virtual ActiveType  theta()                   = 0;
	    virtual ActiveType  m()                       = 0;
	    virtual ActiveType  eta(    const DateType t) = 0;
	    virtual ActiveType  z0()                      = 0;
	    virtual ActiveType  rho()                     = 0;
		// helper functions for vol averaging

/*
Function A_CIR(c1, c2, dt, z0, theta, eta)
    gamma = (theta ^ 2 + 2 * eta ^ 2 * c2) ^ 0.5
    t1 = theta * z0 / eta ^ 2 * (theta + gamma) * dt
    t2 = 1 + (theta + gamma + c1 * eta ^ 2) * (Exp(gamma * dt) - 1) / 2 / gamma
    If t2 <= 0 Then
        A_CIR = "Error! Log of negative argument."
        Exit Function
    End If
    A_CIR = t1 - 2 * theta * z0 / eta ^ 2 * Log(t2)
End Function

Function B_CIR(c1, c2, dt, z0, theta, eta)
    gamma = (theta ^ 2 + 2 * eta ^ 2 * c2) ^ 0.5
    numer = (2 * c2 - theta * c1) * (1 - Exp(-gamma * dt)) + gamma * c1 * (1 + Exp(-gamma * dt))
    denum = (theta + gamma + c1 * eta ^ 2) * (1 - Exp(-gamma * dt)) + 2 * gamma * Exp(-gamma * dt)
    B_CIR = numer / denum
End Function
*/

		// averaging formula implementations
        virtual ActiveType  averageLambda( const DateType T) = 0;
	    virtual ActiveType  averageB     ( const DateType T) = 0;
	    virtual ActiveType  averageEta   ( const DateType T) = 0;


		// embeded classes

		class PieceWiseConstant : TemplateTimeDependentStochVolModel {
		protected:
		    // check for dimensions
		    bool isConsistent_;
		    // time grid
		    std::vector<DateType>     times_;
            // time-dependent model parameters
            std::vector<ActiveType>   lambda_;
	        std::vector<ActiveType>   b_;
	        std::vector<ActiveType>   eta_;
            // time-homogeneous model parameters
		    ActiveType                L_;
	        ActiveType                theta_;
	        ActiveType                m_;
	        ActiveType                z0_;
	        ActiveType                rho_;
			// helper
			inline size_t idx( const DateType t ) { return TemplateAuxilliaries::idx(times_,t); }
		public:
			// constructor
			PieceWiseConstant( const std::vector<DateType>&    times,
				               const std::vector<ActiveType>&  lambda,
							   const std::vector<ActiveType>&  b,
							   const std::vector<ActiveType>&  eta,
							   const ActiveType                L,
							   const ActiveType                theta,
							   const ActiveType                m,
							   const ActiveType                z0,
							   const ActiveType                rho )
			    : times_(times), lambda_(lambda), b_(b), eta_(eta), L_(L), theta_(theta), m_(m), z0_(z0), rho_(rho) {
                isConsistent_ = true;
				QL_REQUIRE(times_.size()>0,"TemplateTimeDependentStochVolModel::PieceWiseConstant: non-empty times required");
				for (size_t k=0; k<times_.size()-1; ++k) QL_REQUIRE(times_[k]<times_[k+1],"TemplateTimeDependentStochVolModel::PieceWiseConstant: ascending time-grid required");
				QL_REQUIRE(lambda_.size()==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: lambda dimension mismatch");
				QL_REQUIRE(b_.size()     ==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: b dimension mismatch");
				QL_REQUIRE(eta_.size()   ==times_.size(),"TemplateTimeDependentStochVolModel::PieceWiseConstant: eta dimension mismatch");
			}
			// inspectors
			virtual ActiveType  lambda( const DateType t)  { return lambda_[idx(t)]; }
	        virtual ActiveType  b(      const DateType t)  { return b_[idx(t)];      }
	        virtual ActiveType  eta(    const DateType t)  { return eta_[idx(t)];    }
			virtual ActiveType  L()                        { return L_;              }
			virtual ActiveType  theta()                    { return theta_;          }
			virtual ActiveType  m()                        { return m_;              }
			virtual ActiveType  z0()                       { return z0_;             }
			virtual ActiveType  rho()                      { return rho_;            }

		};

		// averaging using Gauss-Lobatto integration
		class GaussLobatto : public TemplateTimeDependentStochVolModel {
		protected:
			// Gauss Lobatto options
			PassiveType absAccuracy_;
			PassiveType relAccuracy_;
			size_t      maxEvaluations_;
			// ode solver step size
			DateType    dt_;
			// abbreviation
#define _integr_ TemplateAuxilliaries::GaussLobatto<ActiveType>(m_->maxEvaluations(),  m_->absAccuracy(), m_->relAccuracy()).integrate
			//eta averaging Piterbarg, Th. 9.3.6
			struct f1 {
				GaussLobatto* m_;
				DateType s_;
				f1(GaussLobatto* m, DateType s) : m_(m), s_(s) {}
				inline ActiveType operator() ( DateType t ) { return m_->lambda(t)*m_->lambda(t) * exp(-m_->theta()*(t-s_)); }
			};
			struct I1 {
				GaussLobatto* m_;
				DateType T_;
				I1( GaussLobatto* m, DateType T ) : m_(m), T_(T) {}
				inline ActiveType operator() ( DateType s ) { return _integr_(f1(m_,s),s,T_); }
			};
			struct f2 {
				GaussLobatto* m_;
				DateType r_;
				DateType T_;
				f2(GaussLobatto* m, DateType r, DateType T ) : m_(m), r_(r), T_(T) {}
				inline ActiveType operator() ( DateType s ) { return m_->lambda(s)*m_->lambda(s) * exp(-2.0*m_->theta()*(s-r_)) * I1(m_,T_)(s); }
			};
			struct I2 {  // rho_T(r)
				GaussLobatto* m_;
				DateType T_;
				I2( GaussLobatto* m, DateType T ) : m_(m), T_(T) {}
				inline ActiveType operator() ( DateType r ) { return _integr_(f2(m_,r,T_),r,T_); }
			};
			struct f3 {
				GaussLobatto* m_;
				DateType T_;
				f3(GaussLobatto* m, DateType T ) : m_(m), T_(T) {}
				inline ActiveType operator() ( DateType t ) { return m_->eta(t)*m_->eta(t) * I2(m_,T_)(t); }
			};
			struct I3 {  // numerator
				GaussLobatto* m_;
				I3( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f3(m_,T),0,T); }
			};
			struct I4 {  // denumerator
				GaussLobatto* m_;
				I4( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(I2(m_,T),0,T); }
			};
			// b averaging, Piterbarg, (9.36) rearranged for exp{...}
			struct f4 {
				GaussLobatto* m_;
				f4(GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType s ) { return m_->lambda(s) * m_->lambda(s); }
			};
			struct I5 {
				GaussLobatto* m_;
				I5( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType t ) { return _integr_(f4(m_),0,t); }
			};
			struct f5 {
				GaussLobatto* m_;
				DateType      t_, s_;
				f5(GaussLobatto* m, DateType t, DateType s ) : m_(m), t_(t), s_(s) {}
				inline ActiveType operator() ( DateType u ) { return m_->eta(u) * m_->eta(u) * exp(- m_->theta() * (t_-u + s_-u) ); }
			};
			struct I6 {
				GaussLobatto* m_;
				DateType      t_;
				I6( GaussLobatto* m, DateType t ) : m_(m), t_(t) {}
				inline ActiveType operator() ( DateType s ) { return _integr_(f5(m_,t_,s),0,s); }
			};
			struct f6 {
				GaussLobatto* m_;
				DateType      t_;
				f6(GaussLobatto* m, DateType t ) : m_(m), t_(t) {}
				inline ActiveType operator() ( DateType s ) { return m_->lambda(s) * m_->lambda(s) * I6(m_,t_)(s) ; }
			};
			struct I7 {
				GaussLobatto* m_;
				I7( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType t ) { return _integr_(f6(m_,t),0,t); }
			};
			struct f7 { // \hat{v}(t)^2
				GaussLobatto* m_;
				f7(GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType t ) { return m_->z0() * m_->z0() * I5(m_)(t) + m_->z0() * I7(m_)(t); }
			};
			struct f8 { // numerator integrand
				GaussLobatto* m_;
				f8(GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType t ) { return f7(m_)(t) * m_->lambda(t) * m_->lambda(t) * m_->b(t); }
			};
			struct f9 { // denumerator integrand
				GaussLobatto* m_;
				f9(GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType t ) { return f7(m_)(t) * m_->lambda(t) * m_->lambda(t); }
			};
			struct I8 {
				GaussLobatto* m_;
				I8( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f8(m_),0,T); }
			};
			struct I9 {
				GaussLobatto* m_;
				I9( GaussLobatto* m ) : m_(m) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f9(m_),0,T); }
			};
#undef _integr_
			// Riccati ODE for vol averaging
			struct riccati_f {
				GaussLobatto* m_;
				DateType      u_;
				DateType      v_;
				DateType      T_;
				ActiveType    b_;
				riccati_f ( GaussLobatto* m, DateType u, DateType v, DateType T, ActiveType b )
					: m_(m), u_(u), v_(v), T_(T), b_(m->averageB(T)) {}
				inline void operator() ( DateType t, const std::vector<ActiveType>& y, std::vector<ActiveType>& f ) {
				    // assume y.size()==f.size()==2
				    f[0] = - m_->theta() * m_->z0() * y[1];
				    f[1] = (m_->theta() - m_->rho() * m_->eta(t) * b_ * u_ * m_->lambda(t) ) * y[1]
					       - m_->eta(t) * m_->eta(t) * y[1] * y[1] / 2.0
						   - v_ * m_->lambda(t) * m_->lambda(t) ;
				}
			};
		public:
			// constructor
			GaussLobatto( const PassiveType absAccuracy = 1,
			              const PassiveType relAccuracy = 1.0e-4,
			              const size_t      maxEvaluations = 1000,
						  const DateType    dt = 1)
						  : absAccuracy_(absAccuracy), relAccuracy_(relAccuracy), maxEvaluations_(maxEvaluations), dt_(dt) {}
			// inspectors
			inline PassiveType absAccuracy() const { return absAccuracy_; }
			inline PassiveType relAccuracy() const { return relAccuracy_; }
			inline size_t      maxEvaluations() const { return maxEvaluations_; } 
		    // averaging...
			virtual ActiveType  averageEta   ( const DateType T) {
				ActiveType num = I3(this)(T);
				ActiveType den = I4(this)(T);
				return sqrt(num / den);
			}
			virtual ActiveType  averageB     ( const DateType T) {
				ActiveType num = I8(this)(T);
				ActiveType den = I9(this)(T);
				return num / den;
			}
			virtual ActiveType  averageLambda( const DateType T) {
				return 0; // todo...
			}
		};

	    // piecewise constant parameters and numerical integration
        class PWCNumerical : public PieceWiseConstant, GaussLobatto {
		public:
			PWCNumerical( const std::vector<DateType>&    times,
				          const std::vector<ActiveType>&  lambda,
					      const std::vector<ActiveType>&  b,
					      const std::vector<ActiveType>&  eta,
					      const ActiveType                L,
					      const ActiveType                theta,
					      const ActiveType                m,
					      const ActiveType                z0,
				          const ActiveType                rho,
                          const PassiveType               absAccuracy = 1,
			              const PassiveType               relAccuracy = 1.0e-4,
			              const size_t                    maxEvaluations = 1000,
					      const DateType                  dt = 1 )
					      : PieceWiseConstant(times,lambda,b,eta,L,theta,m,z0,rho),
					        GaussLobatto(absAccuracy,relAccuracy,maxEvaluations,dt) {}
						
			// inspectors
			virtual ActiveType  lambda( const DateType t)  { return PieceWiseConstant::lambda(t);    }
	        virtual ActiveType  b(      const DateType t)  { return PieceWiseConstant::b(t) ;        }
	        virtual ActiveType  eta(    const DateType t)  { return PieceWiseConstant::eta(t);       }
			virtual ActiveType  L()                        { return PieceWiseConstant::L();          }
			virtual ActiveType  theta()                    { return PieceWiseConstant::theta();      }
			virtual ActiveType  m()                        { return PieceWiseConstant::m();          }
			virtual ActiveType  z0()                       { return PieceWiseConstant::z0();         }
			virtual ActiveType  rho()                      { return PieceWiseConstant::rho();        }
			// averaging formula implementations
			virtual ActiveType  averageLambda( const DateType T) { return GaussLobatto::averageLambda( T ); }
			virtual ActiveType  averageB     ( const DateType T) { return GaussLobatto::averageB( T );      }
			virtual ActiveType  averageEta   ( const DateType T) { return GaussLobatto::averageEta( T );    }
		};
		  
	};



}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_hestonmodels_hpp */
