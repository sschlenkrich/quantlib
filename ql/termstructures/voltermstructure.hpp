/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Ferdinando Ametrano

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file voltermstructure.hpp
    \brief Volatility term structure
*/

#ifndef quantlib_vol_term_structure_hpp
#define quantlib_vol_term_structure_hpp

#include <ql/termstructure.hpp>
#include <ql/termstructures/volatility/volatilitytype.hpp>


namespace QuantLib {

    //! Volatility term structure
    /*! This abstract class defines the interface of concrete
        volatility structures which will be derived from this one.

    */
    class VolatilityTermStructure : public TermStructure {
      public:
        /*! \name Constructors
            See the TermStructure documentation for issues regarding
            constructors.
        */
        //@{
        /*! \warning term structures initialized by means of this
                     constructor must manage their own reference date
                     by overriding the referenceDate() method.
        */
        VolatilityTermStructure(BusinessDayConvention bdc,
                                const DayCounter& dc = DayCounter(),
								const VolatilityType& volatilityType = Lognormal);
        //! initialize with a fixed reference date
        VolatilityTermStructure(const Date& referenceDate,
                                const Calendar& cal,
                                BusinessDayConvention bdc,
                                const DayCounter& dc = DayCounter(),
								const VolatilityType& volatilityType = Lognormal);
        //! calculate the reference date based on the global evaluation date
        VolatilityTermStructure(Natural settlementDays,
                                const Calendar& cal,
                                BusinessDayConvention bdc,
                                const DayCounter& dc = DayCounter(),
								const VolatilityType& volatilityType = Lognormal);
        //@}
        //! the business day convention used in tenor to date conversion
		virtual BusinessDayConvention businessDayConvention() const;
		virtual VolatilityType volatilityType() const { return volatilityType_; }
        //! period/date conversion
        Date optionDateFromTenor(const Period&) const;
        //! the minimum strike for which the term structure can return vols
        virtual Rate minStrike() const = 0;
        //! the maximum strike for which the term structure can return vols
        virtual Rate maxStrike() const = 0;
      protected:
        //! strike-range check
        void checkStrike(Rate strike,
                         bool extrapolate) const;
      private:
        BusinessDayConvention bdc_;
        VolatilityType volatilityType_;
    };

    // inline definitions

    inline BusinessDayConvention
    VolatilityTermStructure::businessDayConvention() const {
        return bdc_;
    }

    inline Date
    VolatilityTermStructure::optionDateFromTenor(const Period& p) const {
        // swaption style
        return calendar().advance(referenceDate(),
                                  p,
                                  businessDayConvention());
    }
}

#endif
