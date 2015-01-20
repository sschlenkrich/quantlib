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
#include <ql/experimental/template/auxilliaries/solver1d.hpp>
#include <ql/experimental/template/templatestochasticprocess.hpp>



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
        inline ActiveType kappa() const { return kappa_; }
        inline ActiveType theta() const { return theta_; }
        inline ActiveType sigma() const { return sigma_; }
        inline ActiveType rho()   const { return rho_;   }
        inline ActiveType v0()    const { return v0_;    }
        // maths
        inline bool fellerConstraint() {
            return (sigma >= 0.0 && sigma*sigma < 2.0*kappa*theta);
        }
        // undiscounted expectation of vanilla payoff
        inline ActiveType vanillaOption(const PassiveType forwardPrice,
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
        inline ActiveType vanillaOption(const PassiveType forwardPrice,
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
	class TemplateTimeDependentStochVolModel : public TemplateStochasticProcess<DateType,PassiveType,ActiveType> {
	public:
		// abstract inspectors
        virtual ActiveType  lambda( const DateType t) = 0;
	    virtual ActiveType  b(      const DateType t) = 0;
	    virtual ActiveType  L()                       = 0;
	    virtual ActiveType  theta()                   = 0;
	    virtual ActiveType  m()                       = 0;
	    virtual ActiveType  eta(    const DateType t) = 0;
	    virtual ActiveType  z0()                      = 0;
	    virtual ActiveType  rho()                     = 0;
	    virtual ActiveType  S0()                      = 0;
		// helper functions for vol averaging, Piterbarg, 10.2.4
		inline static ActiveType A_CIR ( ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt ) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
            ActiveType t1 = theta*z0/eta/eta * (theta + gamma)*dt;
            ActiveType t2 = 1.0 + (theta + gamma + c1 * eta*eta) * (exp(gamma * dt) - 1.0) / 2.0 / gamma;
            QL_REQUIRE(t2>0,"A_CIR: require positive log()-argument");
            return t1 - 2.0*theta*z0/eta/eta*log(t2);
		}

		inline static ActiveType B_CIR ( ActiveType c1, ActiveType c2, ActiveType z0, ActiveType theta, ActiveType eta, DateType dt ) {
            ActiveType gamma = sqrt((theta*theta + 2.0 * eta*eta * c2));
			ActiveType emGdt = exp(-gamma * dt);
            ActiveType numer = (2.0*c2 - theta*c1)*(1.0 - emGdt) + gamma*c1*(1.0 + emGdt);
            ActiveType denum = (theta + gamma + c1*eta*eta) * (1.0 - emGdt) + 2.0*gamma*emGdt;
			return numer / denum;
		}		

		// abstract averaging formula definitions
        inline virtual ActiveType  averageLambda( const DateType T) = 0;
	    inline virtual ActiveType  averageB     ( const DateType T) = 0;
	    inline virtual ActiveType  averageEta   ( const DateType T) = 0;

        // undiscounted expectation of vanilla payoff
        inline ActiveType vanillaOption(const PassiveType forwardPrice,
                                        const PassiveType strikePrice,
                                        const DateType    term,
                                        const int         callOrPut,
                                        const PassiveType accuracy,
                                        const size_t      maxEvaluations) {
				TemplateStochVolModel<DateType,PassiveType,ActiveType> model(averageLambda(term), averageB(term),L(),theta(),m(),averageEta(term),z0(),rho());
			    return model.vanillaOption( forwardPrice, strikePrice, term, callOrPut, accuracy, maxEvaluations );
		}

		// stochastic process interface
		// dimension of X
		inline virtual size_t size() { return 2; }
		// stochastic factors of x and z (maybe distinguish if trivially eta=0)
		inline virtual size_t factors() { return 2; }
		// initial values for simulation
		inline virtual VecP initialValues() {
			VecP X(2);
			X[0] = S0();
			X[1] = z0();
			return X;
		}
		// a[t,X(t)]
		inline virtual VecA drift( const DateType t, const VecA& X) {
			VecA a(2);
			// S-variable drift-less
			a[0] = 0.0;
            // z-variable theta [ m - z(t)^+ ]  (full truncation)
			a[1] = theta()*(m() - ((X[1]>0)?(X[1]):(0.0)));
			return a;
		}
		// b[t,X(t)]
		inline virtual MatA diffusion( const DateType t, const VecA& X) {
			MatA B(2);
			B[0].resize(2);
			B[1].resize(2);
			ActiveType sqrtz = ( (X[1]>0) ? (sqrt(X[1])) : (0.0) );   // full truncation
			// S-variable lambda(t) [ b(t) S(t) + (1-b(t)) L ] sqrt[z(t)] dW(t)
			B[0][0] = lambda(t) * (b(t) * X[0] + (1.0-b(t)) * L()) * sqrtz;
			B[0][1] = 0.0;
			// z-variable
			B[1][0] = sqrt(rho())*eta(t)*sqrtz;
			B[1][1] = sqrt(1-rho())*eta(t)*sqrtz;
			// finished
			return B;
		}


		// embeded classes

		class PieceWiseConstant  {
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
	        ActiveType                S0_;
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
							   const ActiveType                rho,
							   const ActiveType                S0 )
			    : times_(times), lambda_(lambda), b_(b), eta_(eta), L_(L), theta_(theta), m_(m), z0_(z0), rho_(rho), S0_(S0) {
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
			virtual ActiveType  S0()                       { return S0_;            }

		};

		// averaging assuming piecewise constant parameters on ( T_k-1, T_k ]
		class MidPointIntegration {
		protected:
			// reference to model
			TemplateTimeDependentStochVolModel* model_;
			// time grid for integration
			std::vector<DateType>               times_;
			inline std::vector<DateType> getTimes( const DateType T ) {
				std::vector<DateType> times;
				times.push_back(0.0);
				for (size_t k=0; k<times_.size(); ++k) {
					if (times_[k]>0.0 && times_[k]<T) times.push_back(times_[k]);
					if (times_[k]>=T) break;
				}
				times.push_back(T);
				return times;
			}

			class AverageLambdaObjective {
			protected:
				ActiveType z0_;
				ActiveType theta_;
				ActiveType eta_;
				DateType   dt_;
				ActiveType target_;
			public:
				AverageLambdaObjective ( const ActiveType z0,
				                         const ActiveType theta,
				                         const ActiveType eta,
				                         const DateType   dt,
				                         const ActiveType target )
										 : z0_(z0), theta_(theta), eta_(eta), dt_(dt), target_(target) {}
				ActiveType operator() ( const ActiveType avLambda2c ) {
					ActiveType A = A_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
					ActiveType B = B_CIR(0, -avLambda2c, z0_, theta_, eta_, dt_);
					ActiveType res = A - B*z0_ - target_;
					return res;
				}
			};


		public:
			// constructor
			MidPointIntegration ( TemplateTimeDependentStochVolModel* model,
				                  std::vector<DateType>               times )
								  : model_(model), times_(times) { }

		    // averaging formula implementations
	        virtual ActiveType  averageEta   ( const DateType T) {
				std::vector<DateType> times(getTimes(T));
				std::vector<ActiveType> f(times.size());
				std::vector<ActiveType> w(times.size()-1);
				f[f.size()-1] = 0.0;
				for (size_t k=f.size()-1; k>0; --k) {
					ActiveType lambda = model_->lambda((times[k-1]+times[k])/2.0);
					ActiveType tmp    = exp(-model_->theta()*(times[k]-times[k-1]));
					f[k-1] = lambda*lambda/model_->theta()*(1.0 - tmp) + tmp*f[k];
				}
				ActiveType sum = 0.0;
				for (size_t k=0; k<w.size(); ++k) {
					ActiveType lambda  = model_->lambda((times[k]+times[k+1])/2.0);
					ActiveType lambda2 = lambda*lambda;
					ActiveType theta   = model_->theta();
					ActiveType theta2  = theta*theta;
					w[k] = (f[k+1]*f[k+1] - f[k]*f[k])/2.0/theta +
						   (f[k+1] - f[k])*lambda2/theta2 +
						   (times[k+1]-times[k])*lambda2*lambda2/theta2;
					w[k] *= 0.5;
					sum  += w[k];
				}
				ActiveType eta2=0.0;
				for (size_t k=0; k<w.size(); ++k) {
					ActiveType eta = model_->eta((times[k]+times[k+1])/2.0);
					eta2 += w[k] * eta * eta;
				}
				eta2 = eta2 / sum;
				return sqrt(eta2);
			}

            virtual ActiveType  averageB( const DateType T) {
				std::vector<DateType> times(getTimes(T));
				std::vector<ActiveType> w(times.size()-1);
				ActiveType theta = model_->theta();
				ActiveType z0 = model_->z0();
				ActiveType sumLambda2dT = 0.0;
				ActiveType S1=0.0, S2=0.0, S3=0.0;
				ActiveType sum = 0.0;
				for (size_t k=0; k<w.size(); ++k) {
					// v1
					ActiveType lambda    = model_->lambda((times[k]+times[k+1])/2.0);
					ActiveType lambda2   = lambda*lambda;
					ActiveType eta       = model_->eta((times[k]+times[k+1])/2.0);
					ActiveType lambda2dT = lambda2*(times[k+1]-times[k]);
					ActiveType v1        = z0*z0*(times[k+1]-times[k])*(lambda2dT/2.0 + sumLambda2dT);
					sumLambda2dT        += lambda2dT;
					// v3, v4, v5
					ActiveType expmThdT  = exp(-theta*(times[k+1]-times[k]));
					ActiveType v3        = z0/theta*(S2+S3)*(1.0-expmThdT);
					ActiveType v4        = (times[k+1]-times[k]) - (1.0-expmThdT)/theta -
						                   (1.0-expmThdT)*(1.0-expmThdT)/2.0/theta;
					ActiveType laEtaTh   = lambda*eta/theta;
					v4                  *= z0*laEtaTh*laEtaTh/2.0;
					ActiveType theta2    = theta*theta;
					ActiveType v5        = z0*lambda2/theta2/2.0*S1*(1.0-expmThdT)*(1.0-expmThdT);
					// updating S1, S2, S3
					S3                   = expmThdT * ( S3 + S1*lambda2/theta*(1.0-expmThdT) );
					S2                   = expmThdT*S2 + laEtaTh*laEtaTh/2.0*(1.0-expmThdT)*(1.0-expmThdT);
					S1                   = expmThdT*expmThdT*S1 + eta*eta/theta/2.0*(1.0 - expmThdT*expmThdT);
					// gathering things together...
					w[k]                 = lambda2 * ( v1 + v3 + v4 + v5 );
					sum                 += w[k];
				}
				ActiveType b=0;
				for (size_t k=0; k<w.size(); ++k) {
					b += w[k] * model_->b((times[k]+times[k+1])/2.0);
				}
				b = b / sum;
				return b;
			}
			
			virtual ActiveType  averageLambda     ( const DateType T) {
				ActiveType b   = averageB(T);
				ActiveType eta = averageEta(T);  // maybe better use time-dep eta
				std::vector<DateType> times(getTimes(T));
				// c = h''(zeta) / h'(zeta)
				ActiveType zeta = 0.0;
				for (size_t k=0; k<times.size()-1; ++k) {
					ActiveType lambda    = model_->lambda((times[k]+times[k+1])/2.0);
					zeta += lambda*lambda*(times[k+1]-times[k]);
				}
				ActiveType avLambda2 = zeta / (times[times.size()-1] - times[0]);
				zeta *= model_->z0();
				ActiveType c = -(b*b/4.0 + 1.0/zeta)/2.0;
				// Psi_{z lambda^2}
				ActiveType A = 0.0, B = 0.0;
				for (size_t k=0; k<times.size()-1; ++k) {
					DateType  t = (times[k]+times[k+1])/2.0;
					DateType dt = (times[k+1]-times[k]);
					ActiveType lambda = model_->lambda(t);
					A = A + A_CIR(B, -c*lambda*lambda, model_->z0(), model_->theta(), eta, dt);
					B = B_CIR(B, -c*lambda*lambda, model_->z0(), model_->theta(), eta, dt);
				}
				ActiveType target = A - B*model_->z0();
				AverageLambdaObjective f(model_->z0(), model_->theta(), eta, T, target);
				ActiveType avLambda2c = avLambda2 * c;
				avLambda2c = TemplateAuxilliaries::solve1d<ActiveType>(f, 1.0e-8, avLambda2c, avLambda2c, 10);
				ActiveType avLambda = sqrt(avLambda2c / c);
				return avLambda;
			};
		};

		// averaging using Gauss-Lobatto integration
		class GaussLobatto  {
		protected:
			// reference to model
			TemplateTimeDependentStochVolModel* model_;
			// Gauss Lobatto options
			PassiveType absAccuracy_;
			PassiveType relAccuracy_;
			size_t      maxEvaluations_;
			// ode solver step size
			DateType    dt_;
			// abbreviation
#define _integr_ TemplateAuxilliaries::GaussLobatto<ActiveType>(g_->maxEvaluations(),  g_->absAccuracy(), g_->relAccuracy()).integrate
			//eta averaging Piterbarg, Th. 9.3.6
			struct f1 {
				GaussLobatto* g_;
				DateType s_;
				f1(GaussLobatto* g, DateType s) : g_(g), s_(s) {}
				inline ActiveType operator() ( DateType t ) { return g_->model()->lambda(t)*g_->model()->lambda(t) * exp(-g_->model()->theta()*(t-s_)); }
			};
			struct I1 {
				GaussLobatto* g_;
				DateType T_;
				I1( GaussLobatto* g, DateType T ) : g_(g), T_(T) {}
				inline ActiveType operator() ( DateType s ) { return _integr_(f1(g_,s),s,T_); }
			};
			struct f2 {
				GaussLobatto* g_;
				DateType r_;
				DateType T_;
				f2(GaussLobatto* g, DateType r, DateType T ) : g_(g), r_(r), T_(T) {}
				inline ActiveType operator() ( DateType s ) { return g_->model()->lambda(s)*g_->model()->lambda(s) * exp(-2.0*g_->model()->theta()*(s-r_)) * I1(g_,T_)(s); }
			};
			struct I2 {  // rho_T(r)
				GaussLobatto* g_;
				DateType T_;
				I2( GaussLobatto* g, DateType T ) : g_(g), T_(T) {}
				inline ActiveType operator() ( DateType r ) { return _integr_(f2(g_,r,T_),r,T_); }
			};
			struct f3 {
				GaussLobatto* g_;
				DateType T_;
				f3(GaussLobatto* g, DateType T ) : g_(g), T_(T) {}
				inline ActiveType operator() ( DateType t ) { return g_->model()->eta(t)*g_->model()->eta(t) * I2(g_,T_)(t); }
			};
			struct I3 {  // numerator
				GaussLobatto* g_;
				I3( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f3(g_,T),0,T); }
			};
			struct I4 {  // denumerator
				GaussLobatto* g_;
				I4( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(I2(g_,T),0,T); }
			};
			// b averaging, Piterbarg, (9.36) rearranged for exp{...}
			struct f4 {
				GaussLobatto* g_;
				f4(GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType s ) { return g_->model()->lambda(s) * g_->model()->lambda(s); }
			};
			struct I5 {
				GaussLobatto* g_;
				I5( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType t ) { return _integr_(f4(g_),0,t); }
			};
			struct f5 {
				GaussLobatto* g_;
				DateType      t_, s_;
				f5(GaussLobatto* g, DateType t, DateType s ) : g_(g), t_(t), s_(s) {}
				inline ActiveType operator() ( DateType u ) { return g_->model()->eta(u) * g_->model()->eta(u) * exp(- g_->model()->theta() * (t_-u + s_-u) ); }
			};
			struct I6 {
				GaussLobatto* g_;
				DateType      t_;
				I6( GaussLobatto* g, DateType t ) : g_(g), t_(t) {}
				inline ActiveType operator() ( DateType s ) { return _integr_(f5(g_,t_,s),0,s); }
			};
			struct f6 {
				GaussLobatto* g_;
				DateType      t_;
				f6(GaussLobatto* g, DateType t ) : g_(g), t_(t) {}
				inline ActiveType operator() ( DateType s ) { return g_->model()->lambda(s) * g_->model()->lambda(s) * I6(g_,t_)(s) ; }
			};
			struct I7 {
				GaussLobatto* g_;
				I7( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType t ) { return _integr_(f6(g_,t),0,t); }
			};
			struct f7 { // \hat{v}(t)^2
				GaussLobatto* g_;
				f7(GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType t ) { return g_->model()->z0() * g_->model()->z0() * I5(g_)(t) + g_->model()->z0() * I7(g_)(t); }
			};
			struct f8 { // numerator integrand
				GaussLobatto* g_;
				f8(GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType t ) { return f7(g_)(t) * g_->model()->lambda(t) * g_->model()->lambda(t) * g_->model()->b(t); }
			};
			struct f9 { // denumerator integrand
				GaussLobatto* g_;
				f9(GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType t ) { return f7(g_)(t) * g_->model()->lambda(t) * g_->model()->lambda(t); }
			};
			struct I8 {
				GaussLobatto* g_;
				I8( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f8(g_),0,T); }
			};
			struct I9 {
				GaussLobatto* g_;
				I9( GaussLobatto* g ) : g_(g) {}
				inline ActiveType operator() ( DateType T ) { return _integr_(f9(g_),0,T); }
			};
#undef _integr_
			// Riccati ODE for vol averaging
			struct riccati_f {
				GaussLobatto* g_;
				DateType      u_;
				DateType      v_;
				DateType      T_;
				ActiveType    b_;
				riccati_f ( GaussLobatto* g, DateType u, DateType v, DateType T, ActiveType b )
					: g_(g), u_(u), v_(v), T_(T), b_(m->averageB(T)) {}
				inline void operator() ( DateType t, const std::vector<ActiveType>& y, std::vector<ActiveType>& f ) {
				    // assume y.size()==f.size()==2
				    f[0] = - g_->model()->theta() * g_->model()->z0() * y[1];
				    f[1] = (g_->model()->theta() - g_->model()->rho() * g_->model()->eta(t) * b_ * u_ * g_->model()->lambda(t) ) * y[1]
					       - g_->model()->eta(t) * g_->model()->eta(t) * y[1] * y[1] / 2.0
						   - v_ * g_->model()->lambda(t) * g_->model()->lambda(t) ;
				}
			};
		public:
			// constructor
			GaussLobatto( TemplateTimeDependentStochVolModel* model,
				          const PassiveType absAccuracy = 1,
			              const PassiveType relAccuracy = 1.0e-4,
			              const size_t      maxEvaluations = 1000,
						  const DateType    dt = 1)
						  : model_(model), absAccuracy_(absAccuracy), relAccuracy_(relAccuracy), maxEvaluations_(maxEvaluations), dt_(dt) {}
			// inspectors
			inline TemplateTimeDependentStochVolModel* model() { return model_; }
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
        class PWCNumerical : public TemplateTimeDependentStochVolModel {
		private:
			boost::shared_ptr<PieceWiseConstant> pwc_;
			boost::shared_ptr<GaussLobatto>      gl_;
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
						  const ActiveType                S0,
                          const PassiveType               absAccuracy = 1,
			              const PassiveType               relAccuracy = 1.0e-4,
			              const size_t                    maxEvaluations = 1000,
					      const DateType                  dt = 1 )
 		        :    pwc_(new PieceWiseConstant(times,lambda,b,eta,L,theta,m,z0,rho,S0)),
				     gl_(new GaussLobatto(this,absAccuracy,relAccuracy,maxEvaluations,dt)) {}						
			// inspectors
			virtual ActiveType  lambda( const DateType t)  { return pwc_->lambda(t);    }
	        virtual ActiveType  b(      const DateType t)  { return pwc_->b(t) ;        }
	        virtual ActiveType  eta(    const DateType t)  { return pwc_->eta(t);       }
			virtual ActiveType  L()                        { return pwc_->L();          }
			virtual ActiveType  theta()                    { return pwc_->theta();      }
			virtual ActiveType  m()                        { return pwc_->m();          }
			virtual ActiveType  z0()                       { return pwc_->z0();         }
			virtual ActiveType  rho()                      { return pwc_->rho();        }
			virtual ActiveType  S0()                       { return pwc_->S0();         }
			// averaging formula implementations
			virtual ActiveType  averageLambda( const DateType T) { return gl_->averageLambda( T ); }
			virtual ActiveType  averageB     ( const DateType T) { return gl_->averageB( T );      }
			virtual ActiveType  averageEta   ( const DateType T) { return gl_->averageEta( T );    }
		};

			    // piecewise constant parameters and numerical integration
        class PWCAnalytical : public TemplateTimeDependentStochVolModel {
		private:
			boost::shared_ptr<PieceWiseConstant>    pwc_;
			boost::shared_ptr<MidPointIntegration>  mp_;
		public:
			PWCAnalytical(const std::vector<DateType>&    times,
				          const std::vector<ActiveType>&  lambda,
					      const std::vector<ActiveType>&  b,
					      const std::vector<ActiveType>&  eta,
					      const ActiveType                L,
					      const ActiveType                theta,
					      const ActiveType                m,
					      const ActiveType                z0,
				          const ActiveType                rho,
						  const ActiveType                S0,
                          const PassiveType               absAccuracy = 1,
			              const PassiveType               relAccuracy = 1.0e-4,
			              const size_t                    maxEvaluations = 1000,
					      const DateType                  dt = 1 )
 		        :    pwc_(new PieceWiseConstant(times,lambda,b,eta,L,theta,m,z0,rho,S0)),
				     mp_(new MidPointIntegration(this,times)) {}						
			// inspectors
			virtual ActiveType  lambda( const DateType t)  { return pwc_->lambda(t);    }
	        virtual ActiveType  b(      const DateType t)  { return pwc_->b(t) ;        }
	        virtual ActiveType  eta(    const DateType t)  { return pwc_->eta(t);       }
			virtual ActiveType  L()                        { return pwc_->L();          }
			virtual ActiveType  theta()                    { return pwc_->theta();      }
			virtual ActiveType  m()                        { return pwc_->m();          }
			virtual ActiveType  z0()                       { return pwc_->z0();         }
			virtual ActiveType  rho()                      { return pwc_->rho();        }
			virtual ActiveType  S0()                       { return pwc_->S0();         }
			// averaging formula implementations
			virtual ActiveType  averageLambda( const DateType T) { return mp_->averageLambda( T ); }
			virtual ActiveType  averageB     ( const DateType T) { return mp_->averageB( T );      }
			virtual ActiveType  averageEta   ( const DateType T) { return mp_->averageEta( T );    }
		};


	};



}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_hestonmodels_hpp */
