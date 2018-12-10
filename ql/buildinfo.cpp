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

/*! \file buildinfo.cpp
    \brief Incorporate git build information
*/

#include <ql/buildinfo.hpp>

namespace QuantLib {
    
    std::vector<std::string> buildinfo(){
		std::vector<std::string> info;
		info.push_back("Git repo:   " + std::string(GIT_REPO));
		info.push_back("Git branch: " + std::string(GIT_BRANCH));
		info.push_back("Git commit: " + std::string(GIT_COMMIT_ID));
		info.push_back("Git status: " + std::string(GIT_WORKING_DIR_STATUS));
		info.push_back("Build time: " + std::string(GIT_BUILD_INFO_NOW));
		return info;
	}
    
}
