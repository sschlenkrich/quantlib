/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015, Sebastian Schlenkrich

*/

/*! \file commoditymodels.hpp
    \brief bindings of templatised commodity models
*/


#ifndef quantlib_templatecommoditymodels_hpp
#define quantlib_templatecommoditymodels_hpp

#include <ql/experimental/template/commodity/template2fnormalmodel.hpp>
#include <ql/experimental/template/commodity/template2flognormalmodel.hpp>

namespace QuantLib {

	// basic binding of template parameters
	typedef Template2FMeanReversionModel<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FMeanReversionModel;

	// basic binding of template parameters
	typedef Template2FNormalModel<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FNormalModel;

	// basic binding of template parameters
	typedef Template2FLognormalModel<QuantLib::Time,QuantLib::Real,QuantLib::Real> Real2FLognormalModel;

}

#endif  /* ifndef quantlib_templatecommoditymodels_hpp */
