/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file montecarlomodells.hpp
    \brief bindings of templatised montecarlo simulation and pricing
*/


#ifndef quantlib_montecarlomodells_hpp
#define quantlib_montecarlomodells_hpp

//#include <ql/experimental/templatehullwhite/adtageo/adtageo.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>

#include <ql/experimental/template/templatestochasticprocess.hpp>
#include <ql/experimental/template/montecarlo/templatemcsimulation.hpp>
#include <ql/experimental/template/montecarlo/templatemcpayoff.hpp>


namespace QuantLib {

	// basic binding of template parameters

	typedef TemplateMCSimulation<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCSimulation;

	typedef TemplateMCPayoff< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCPayoff;

	typedef RealMCPayoff::Pricer RealMCPayoffPricer;

}

#endif  /* ifndef quantlib_templatequasigaussianmodels_hpp */
