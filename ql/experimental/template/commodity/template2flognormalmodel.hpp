/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/



#ifndef quantlib_template2fnormalmodel_hpp
#define quantlib_template2fnormalmodel_hpp

#include <ql/experimental/template/commodity/template2fmeanreversionmodel.hpp>



#define _MIN_( a, b ) ( (a) < (b) ? (a) : (b) )
#define _MAX_( a, b ) ( (a) > (b) ? (a) : (b) )

namespace QuantLib {

    // 2-factor mean reverting normal model
	//
	//  ln X(t) = phi(t) + Y(t) + Z(t), X(0)=Y(0)=0
	//    dY(t) = -a Y(t) dt  +  sigma(t) dW_Y(t)
	//    dZ(t) = -b Z(t) dt  +    eta(t) dW_Z(t)
	//    dW_Y(t) dW_Z(t) = rho dt
	//
	template <class DateType, class PassiveType, class ActiveType>
	class Template2FLognormalModel : public Template2FMeanReversionModel<DateType,PassiveType,ActiveType> {
	
	public:
		// constructor
		
		Template2FLognormalModel( const Handle<IndexTermStructure>& futureTS,
		                       const VecD&                          times,
							   const VecA&                          sigma,
							   const VecA&                          eta,
							   const PassiveType                     a,
							   const PassiveType                     b,
							   const PassiveType                     rho 
							   )
		: Template2FMeanReversionModel(futureTS,times,sigma,eta,a,b,rho) {
			// check for valid parameter inputs
		}
	
		// analytic formulas in base class

		inline virtual const ActiveType phi(const DateType t) const {
			ActiveType var = varianceY(0,t) + varianceZ(0,t) + 2.0*rho_*covarianceYZ(0,t);
			return log(futureTS_->value(t)) - var/2.0;
		}


		// basic instruments

		// future expectation
        inline virtual ActiveType futureAsset(const DateType t, const DateType T, const ActiveType Y, const ActiveType Z) {
			ActiveType mu  = phi(T) + exp(-a_*(T-t))*Y + exp(-b_*(T-t))*Z;
			ActiveType var = varianceY(t,T) + varianceZ(t,T) + 2.0*rho_*covarianceYZ(t,T);
			return exp(mu + var/2);
		}

        inline virtual ActiveType averageFuture ( const VecD& settlementTimes, const VecP& settlementWeights) {
			ActiveType fut = 0.0;
			size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
			for (size_t k=0; k<N; ++k) fut += settlementWeights[k] * futureAsset(0.0,settlementTimes[k],0.0,0.0);
			return fut;
		}

		// calculate approximate log-variance as input to Black formula below
		// NOTE: approximation quality needs to be verified
        inline virtual ActiveType varianceAverageFuture ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights) {
			ActiveType  var=0.0;
			PassiveType B = 0, C = 0, B2 = 0, C2 = 0, BC = 0;
			PassiveType weight=0.0;
			size_t N = _MIN_(settlementTimes.size(),settlementWeights.size());
			for (size_t k=0; k<N; ++k) {
				B = exp(-a_*(settlementTimes[k]-expiryTime));
				C = exp(-b_*(settlementTimes[k]-expiryTime));
				B2 += settlementWeights[k]*B*B;
				C2 += settlementWeights[k]*C*C;
				BC += settlementWeights[k]*B*C;
				weight += settlementWeights[k];
				//var += settlementWeights[k]*(B*B*varianceY(0.0,expiryTime) + C*C*varianceZ(0.0,expiryTime) + 2.0*rho_*B*C*covarianceYZ(0.0,expiryTime));
			}
			var = B2*varianceY(0.0,expiryTime) + C2*varianceZ(0.0,expiryTime) + 2.0*rho_*BC*covarianceYZ(0.0,expiryTime);
			var /= weight;
			return var;
		}

		inline virtual ActiveType vanillaOption ( const DateType expiryTime, const VecD& settlementTimes, const VecP& settlementWeights, PassiveType strike, int callOrPut) {
			ActiveType fut = averageFuture(settlementTimes,settlementWeights);
			ActiveType var = varianceAverageFuture(expiryTime,settlementTimes,settlementWeights);
			ActiveType pv  = TemplateAuxilliaries::Black76(fut,strike,sqrt(var),1.0,callOrPut);
			return pv;
		}

	};

}

#undef _MIN_
#undef _MAX_

#endif  /* quantlib_template2fnormalmodel_hpp */
