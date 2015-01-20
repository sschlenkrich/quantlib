/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file quasigaussianmodels.hpp
    \brief bindings of templatised quasi Gaussian models
*/


#ifndef quantlib_templatequasigaussianmodels_hpp
#define quantlib_templatequasigaussianmodels_hpp

//#include <ql/experimental/templatehullwhite/adtageo/adtageo.hpp>
//#include <ql/experimental/template/auxilliaries/MinimADVariable2.hpp>
//#include <ql/experimental/template/templatestochasticprocess.hpp>
#include <ql/experimental/template/qgaussian/templatequasigaussian.hpp>
#include <ql/experimental/template/montecarlo/templatemcsimulation.hpp>
#include <ql/experimental/template/montecarlo/templatemcpayoff.hpp>


namespace QuantLib {

	typedef TemplateStochasticProcess<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealStochasticProcess;

	// basic binding of template parameters
	typedef TemplateQuasiGaussianModel<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealQuasiGaussianModel;

	typedef TemplateMCSimulation<QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCSimulation;

	typedef TemplateMCPayoff< QuantLib::Time,QuantLib::Real,QuantLib::Real> RealMCPayoff;

	typedef RealMCPayoff::Pricer RealMCPayoffPricer;


}

#endif  /* ifndef quantlib_templatequasigaussianmodels_hpp */
