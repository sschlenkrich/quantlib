/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2014 Sebastian Schlenkrich

 */

/*! \file basisswapengine.cpp
    \brief Discounting swap engine for basis swaps
*/

#include <ql/experimental/basisswap/basisswapengine.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/utilities/dataformatters.hpp>



namespace QuantLib {

//
	BasisSwapEngine::BasisSwapEngine(     
		const std::vector<Handle<YieldTermStructure>>& discCurves,
        const std::vector<Real>&                       fxForDom,
		boost::optional<bool>                          includeSettlementDateFlows,
        Date                                           settlementDate,
        Date                                           npvDate
		) : discCurves_(discCurves), fxForDom_(fxForDom),
		    includeSettlementDateFlows_(includeSettlementDateFlows),
            settlementDate_(settlementDate), npvDate_(npvDate) {
        // consistency checks are postponed to calculate() since curve handles may be empty at this time
		for (Size i=0; i<discCurves_.size(); ++i) registerWith(discCurves_[i]);
    }

    void BasisSwapEngine::calculate() const {

        // check for consistency of yield curves
        QL_REQUIRE(!discCurves_.empty(),
                   "discounting term structure handle vector is empty");
        Date refDate = discCurves_[0]->referenceDate();
		for (Size i=1; i<discCurves_.size(); ++i)
			QL_REQUIRE(refDate==discCurves_[i]->referenceDate(),
                   "discounting term structure reference dates differ");

        results_.value = 0.0;
        results_.errorEstimate = Null<Real>();

        Date settlementDate = settlementDate_;
        if (settlementDate_==Date()) {
            settlementDate = refDate;
        } else {
            QL_REQUIRE(settlementDate>=refDate,
                       "settlement date (" << settlementDate << ") before "
                       "discount curve reference date (" << refDate << ")");
        }

        results_.valuationDate = npvDate_;
        if (npvDate_==Date()) {
            results_.valuationDate = refDate;
        } else {
            QL_REQUIRE(npvDate_>=refDate,
                       "npv date (" << npvDate_  << ") before "
                       "discount curve reference date (" << refDate << ")");
        }
        results_.npvDateDiscount = discCurves_[0]->discount(results_.valuationDate);

        Size n = arguments_.legs.size();
        results_.legNPV.resize(n);
        results_.legBPS.resize(n);
        results_.startDiscounts.resize(n);
        results_.endDiscounts.resize(n);

        bool includeRefDateFlows =
            includeSettlementDateFlows_ ?
            *includeSettlementDateFlows_ :
            Settings::instance().includeReferenceDateEvents();

        for (Size i=0; i<n; ++i) {
            try {
				const YieldTermStructure& discount_ref = **discCurves_[std::min(i,discCurves_.size()-1)];
                CashFlows::npvbps(arguments_.legs[i],
                                  discount_ref,
                                  includeRefDateFlows,
                                  settlementDate,
                                  results_.valuationDate,
                                  results_.legNPV[i],
                                  results_.legBPS[i]);
                results_.legNPV[i] *= arguments_.payer[i] * fxForDom_[std::min(i,fxForDom_.size()-1)];
                results_.legBPS[i] *= arguments_.payer[i] * fxForDom_[std::min(i,fxForDom_.size()-1)];

                if (!arguments_.legs[i].empty()) {
                    Date d1 = CashFlows::startDate(arguments_.legs[i]);
                    if (d1>=refDate)
                        results_.startDiscounts[i] = discCurves_[std::min(i,discCurves_.size()-1)]->discount(d1);
                    else
                        results_.startDiscounts[i] = Null<DiscountFactor>();

                    Date d2 = CashFlows::maturityDate(arguments_.legs[i]);
                    if (d2>=refDate)
                        results_.endDiscounts[i] = discCurves_[std::min(i,discCurves_.size()-1)]->discount(d2);
                    else
                        results_.endDiscounts[i] = Null<DiscountFactor>();
                } else {
                    results_.startDiscounts[i] = Null<DiscountFactor>();
                    results_.endDiscounts[i] = Null<DiscountFactor>();
                }

            } catch (std::exception &e) {
                QL_FAIL(io::ordinal(i+1) << " leg: " << e.what());
            }
            results_.value += results_.legNPV[i];
        }
    }



}


