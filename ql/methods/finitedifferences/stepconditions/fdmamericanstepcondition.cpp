/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Andreas Gaida
 Copyright (C) 2008, 2009 Ralph Schreyer
 Copyright (C) 2008, 2009 Klaus Spanderen

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

#include <ql/methods/finitedifferences/operators/fdmlinearoplayout.hpp>
#include <ql/methods/finitedifferences/stepconditions/fdmamericanstepcondition.hpp>

namespace QuantLib {

    FdmAmericanStepCondition::FdmAmericanStepCondition(
            const ext::shared_ptr<FdmMesher> & mesher,
            const ext::shared_ptr<FdmInnerValueCalculator> & calculator)
    : mesher_(mesher),
      calculator_(calculator) {
    }

    void FdmAmericanStepCondition::applyTo(Array& a, Time t) const {
        ext::shared_ptr<FdmLinearOpLayout> layout = mesher_->layout();

        QL_REQUIRE(layout->size() == a.size(),
                   "inconsistent array dimensions");

        const FdmLinearOpIterator endIter = layout->end();

        for (FdmLinearOpIterator iter = layout->begin(); iter != endIter;
            ++iter) {
            const Real innerValue = calculator_->innerValue(iter, t);
            if (innerValue > a[iter.index()]) {
                a[iter.index()] = innerValue;
            }
        }
    }
}
