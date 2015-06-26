/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010, Sebastian Schlenkrich

*/



#ifndef quantlib_qgaussiancalibrator_hpp
#define quantlib_qgaussiancalibrator_hpp

#include <ql/instruments/swaption.hpp>

#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>

#include <ql/experimental/template/auxilliaries/templateauxilliaries.hpp>

#include <ql/experimental/template/basismodel/swaptioncfs.hpp>

#include <ql/experimental/template/qgaussian/quasigaussianmodels.hpp>
#include <ql/experimental/template/montecarlo/montecarlomodells.hpp>


#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {


	// calibrate stochastic volatility model to implied normal volatilities
    class QuasiGaussianModelCalibrator { 

		// the Quasi-Gaussian model aimed to be calibrated
		boost::shared_ptr<RealQuasiGaussianModel> model_;

		// we may use a MC simulation to adjust for closed form vs MC bias
		boost::shared_ptr<RealMCSimulation> mcSimulation_;

		// calibration targets

		// grid of reference (ATM) swaptions as calibration instruments, swaptions_[exercise][swapterm]
		std::vector< std::vector< boost::shared_ptr<Swaption> > > swaptions_;

		// grid of reference swaption stoch vol model parameters
		// user supplied or calculated
		// dimensions according to dimensions of swaptions_
		std::vector< std::vector< Real > > lambda_;
		std::vector< std::vector< Real > > b_;
		std::vector< std::vector< Real > > eta_;

		// constraints for input model parameters
		Real  lambdaMin_, lambdaMax_, bMin_, bMax_, etaMin_, etaMax_;
		// transformation (-inf, +inf) -> (a, b)
        static const Real direct(const Real x, const Real a, const Real b) {
			// return (b-a)*(atan(x)/M_PI + 0.5) + a;
			return TemplateAuxilliaries::direct(x, a, b);
		}
		// transformation (a, b) -> (-inf, +inf)
		static const Real inverse(const Real y, const Real a, const Real b) {
			// return tan( ((y-a)/(b-a)-0.5) * M_PI );
			return TemplateAuxilliaries::inverse(y, a, b);
		}

		class Objective : CostFunction {
		private:
			// reference to access model, swaptions and targets
			QuasiGaussianModelCalibrator  *calibrator_;
			// we work on a copy of the initial model (better safe then sorry)
			boost::shared_ptr<RealQuasiGaussianModel> model_;
			// specify parameters used for optimisation; dimension [modeltimes] x [ d (lambda) + d (b) + 1 (eta) ]
			std::vector< std::vector< bool > > isInput_;
			// specify targets used for optimisation; dimension [exercises] x [ swapterm (lambda) + swapterm (b) + swapterm (eta) ]
			std::vector< std::vector< bool > > isOutput_;
			// count inputs and outputs
			Size inputSize_;
			Size outputSize_;

			// container for single curve swaption cash flow representation and calibration targets; initialised at inception
			class CalibSwaption : public SwaptionCashFlows {
                // cash flows are inherited from base class
				// targets
				Real lambda_, b_, eta_;
			public:
				CalibSwaption ( Real lambda,
					            Real b,
								Real eta,
					            const boost::shared_ptr<Swaption>& swaption,
			                    const Handle<YieldTermStructure>& discountCurve,
						        bool                              contTenorSpread = true )
								: SwaptionCashFlows(swaption,discountCurve,contTenorSpread), lambda_(lambda), b_(b), eta_(eta) {}
				// inspectors
				inline Real lambda() { return lambda_; }
				inline Real b()      { return b_;      }
				inline Real eta()    { return eta_;    }
			};
			std::vector< CalibSwaption > swaptions_;
			
		public:
			const Size inputSize()  const { return inputSize_; }
			const Size outputSize() const { return outputSize_; }

			Objective( QuasiGaussianModelCalibrator              *calibrator,
			           const std::vector< std::vector< bool > >&  isInput,
			           std::vector< std::vector< bool > >&        isOutput );

		    // initialize state X with model parameters and apply inverse transformation
		    Array initialise();
			    // allocate X with inputSize
			    // for all isInput[i]
			    //     for j=0.. d-1, isInput[i][j] ? lambda -> inverseLambda -> X[idx], ++idx
			    //     for j=d..2d-1, isInput[i][j] ? b      -> inverseB      -> X[idx], ++idx
			    //     for j=2d,      isInput[i][j] ? eta    -> inverseEta    -> X[idx], ++idx

			// apply direct transformation and update model parameters
			void update(const Array& X);
			    // for all isInput[i]
			    //     for j=0.. d-1, isInput[i][j] ? X[idx] -> inverseLambda -> lambda, ++idx
			    //     for j=d..2d-1, isInput[i][j] ? X[idx] -> inverseB      -> b     , ++idx
			    //     for j=2d,      isInput[i][j] ? X[idx] -> inverseEta    -> eta   , ++idx
			    //  model->update(lambda,b,eta)

			// CostFunction interface
            virtual Disposable<Array> values(const Array& x) const;
			    // update(x)
			    // allocate Y with outputSize()
			    // for k=0..swaptions_.size() construct RealQGSwaptionModel swpModel[k]
			    // for all isOutput[i], j=0..swapterm-1
			    //     isOutput[i][0*swapterm+j] ? swpModel[k] -> avLambda -> Y[idx], Y[idx] -= lambda[k], ++idx
			    //     isOutput[i][1*swapterm+j] ? swpModel[k] -> avB      -> Y[idx], Y[idx] -= b[k],      ++idx 
			    //     isOutput[i][2*swapterm+j] ? swpModel[k] -> avEta    -> Y[idx], Y[idx] -= eta[k],    ++idx
			    //     ++k

			// 0.5 ||Y||^2
            virtual Real value(const Array& x) const;

		};

	public:

    };

}

#undef _MIN_
#undef _MAX_

#endif  /* ifndef quantlib_qgaussiancalibrator_hpp */
