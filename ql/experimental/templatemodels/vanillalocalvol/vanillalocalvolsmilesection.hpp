/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2017, Sebastian Schlenkrich

*/



#ifndef quantlib_vanillalocalvolsmilesection_hpp
#define quantlib_vanillalocalvolsmilesection_hpp

#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/math/optimization/method.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/termstructures/volatility/smilesection.hpp>

#include <ql/experimental/templatemodels/vanillalocalvol/vanillalocalvolmodels.hpp>

namespace QuantLib {

	class VanillaLocalVolModelSmileSection : public SmileSection {
	private:
		boost::shared_ptr<VanillaLocalVolModel>  model_;
	protected:
		virtual Volatility volatilityImpl(Rate strike) const;
	public:
		VanillaLocalVolModelSmileSection(
			const boost::shared_ptr<VanillaLocalVolModel>&    model,
			const DayCounter&                                 dc = DayCounter(),
			const VolatilityType                              type = Normal,
			const Rate                                        shift = 0.0)
			: model_(model), SmileSection(model->timeToExpiry(), dc, type, shift) {}

		VanillaLocalVolModelSmileSection(
			const Date&                                       expiryDate,
			const Rate&                                       forward,
			const std::vector<Rate>&                          relativeStrikes,
			const std::vector<Volatility>&                    smileVolatilities,
			const Real                                        extrapolationRelativeStrike,
			const Real                                        extrapolationSlope,
			bool                                              vegaWeighted = true,
			const boost::shared_ptr<EndCriteria>&             endCriteria = boost::shared_ptr<EndCriteria>(),
			const boost::shared_ptr<OptimizationMethod>&      method = boost::shared_ptr<OptimizationMethod>(),
			const DayCounter&                                 dc = Actual365Fixed(),
			const Date&                                       referenceDate = Date(),
			const VolatilityType                              type = Normal,
			const Rate                                        shift = 0.0,
			const boost::shared_ptr<VanillaLocalVolModel>&    model = 0,
			const Real                                        minSlope = -1.0,   //  lower boundary for m in calibration
			const Real                                        maxSlope = 1.0,    //  upper boundary for m in calibration
			const Real                                        alpha = 0.0);      //  Tikhonov alpha

		// SmileSection interface 
		virtual Real minStrike() const { return model_->underlyingS().front(); }
		virtual Real maxStrike() const { return model_->underlyingS().back();  }
		virtual Real atmLevel()  const { return model_->forward();             }

		// overload optionPrice() as a basis for implied volatility calculation
		virtual Real optionPrice(Rate strike, Option::Type type = Option::Call, Real discount = 1.0) const;

		// inspector
		inline const boost::shared_ptr<VanillaLocalVolModel>&  model() const { return model_; }
	};
}

#endif  /* ifndef quantlib_vanillalocalvolsmilesection_hpp */
