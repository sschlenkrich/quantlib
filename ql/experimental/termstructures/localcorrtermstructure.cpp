/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C)  2017 Cord Harms

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

#include <ql/experimental/termstructures/localcorrtermstructure.hpp>

namespace QuantLib {

    LocalCorrTermStructure::LocalCorrTermStructure(const std::vector<boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>>& processes, 
												   const boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess>&								   processToCal)
    : CorrelationTermStructureStrike(processes[0]->blackVolatility()->referenceDate(), 
		processes[0]->blackVolatility()->calendar(), 
		processes[0]->blackVolatility()->businessDayConvention(), 
		processes[0]->blackVolatility()->dayCounter()), processes_(processes), processToCal_(processToCal){}

    void LocalCorrTermStructure::localCorr(RealStochasticProcess::MatA& corrMatrix, 
											   const Date& d,
											   const RealStochasticProcess::VecA& X0,
                                               bool extrapolate) {
        
		for (size_t i = 0; i < X0.size(); i++)
		{
			checkRange(d, extrapolate);
			checkStrike(X0[i], i, extrapolate);
		} 
        Time t = timeFromReference(d);
        
		localCorrImpl(corrMatrix, t, X0,extrapolate);
    }

    void LocalCorrTermStructure::localCorr(RealStochasticProcess::MatA& corrMatrix, 
											   Time t,
											   const RealStochasticProcess::VecA& X0,
                                               bool extrapolate) {
		checkRange(t, extrapolate);
		for (size_t i = 0; i < X0.size(); i++)
		{
			checkStrike(X0[i], i, extrapolate);
		}
        localCorrImpl(corrMatrix, t, X0, extrapolate);

		//cap to 1 and floor to -1:
		for (size_t i = 0; i < corrMatrix.size(); i++)
		{
			for (size_t j = i+1; j < corrMatrix.size(); j++)
			{
				if (corrMatrix[i][j] > 1) {
					corrMatrix[i][j] = 1;
					corrMatrix[j][i] = 1;
				}
				if (corrMatrix[i][j] < -1) {
					corrMatrix[i][j] = -1;
					corrMatrix[j][i] = -1;
				}
			}
		}
    }

    void LocalCorrTermStructure::accept(AcyclicVisitor& v) {
        Visitor<LocalCorrTermStructure>* v1 =
            dynamic_cast<Visitor<LocalCorrTermStructure>*>(&v);
        if (v1 != 0)
            v1->visit(*this);
        else
            QL_FAIL("not a local-Correlation term structure visitor");
    }

	const Date& LocalCorrTermStructure::referenceDate() const {
		return processes_[0]->blackVolatility()->referenceDate();
	}

	DayCounter LocalCorrTermStructure::dayCounter() const {
		return processes_[0]->blackVolatility()->dayCounter();
	}

	Date LocalCorrTermStructure::maxDate() const {
		Date minMaxDate = processes_[0]->blackVolatility()->maxDate();
		for (size_t i = 0; i < processes_.size(); i++)
		{
			if (minMaxDate < processes_[i]->blackVolatility()->maxDate()) minMaxDate = processes_[i]->blackVolatility()->maxDate();
		}
		return minMaxDate;
	}

	Real LocalCorrTermStructure::minStrike(Natural ulId) const {
		return processes_[ulId]->blackVolatility()->minStrike();
	}

	Real LocalCorrTermStructure::maxStrike(Natural ulId) const {
		return processes_[ulId]->blackVolatility()->maxStrike();
	}
}
