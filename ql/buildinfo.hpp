/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2018 Sebastian Schlenkrich

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

/*! \file buildinfo.hpp
    \brief Incorporate git build information
*/

#ifndef quantlib_buildinfo_hpp
#define quantlib_buildinfo_hpp

#include <vector>
#include <string>

// use this flag to control external build info use
#define QL_USE_BUILD_INFO
#ifdef QL_USE_BUILD_INFO
    #include <buildinfo.hpp>  // this file is auto-generated
#endif

// we set defaults such that we can use defines

#ifndef GIT_REPO
    #define GIT_REPO "Unknown"
#endif

#ifndef GIT_BRANCH
    #define GIT_BRANCH "Unknown"
#endif

#ifndef GIT_COMMIT_ID
    #define GIT_COMMIT_ID "Unknown"
#endif

#ifndef GIT_COMMIT_ID_SHORT
    #define GIT_COMMIT_ID_SHORT "Unknown"
#endif

#ifndef GIT_WORKING_DIR_STATUS
    #define GIT_WORKING_DIR_STATUS "Unknown\n"
#endif

#ifndef GIT_BUILD_INFO_NOW 
    #define GIT_BUILD_INFO_NOW  "Unknown\n"
#endif


namespace QuantLib {
    
    std::vector<std::string> buildinfo(); 
    
}

#endif  // quantlib_buildinfo_hpp