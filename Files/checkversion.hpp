//
//  check_version.h
//  Master Thesis
//
//  Created by Magnus Mencke on 30/03/2021.
//  Copyright © 2021 Magnus Mencke. All rights reserved.
//

#ifndef check_version_hpp
#define check_version_hpp


#include <iostream>
#include <boost/version.hpp>
#include <ql/version.hpp>


inline void checkVersion() {
    std::cout << "Using Boost version "
    << BOOST_VERSION / 100000     << "."  // major version
    << BOOST_VERSION / 100 % 1000 << "."  // minor version
    << BOOST_VERSION % 100                // patch level
    << std::endl;

    std::cout << "Using QuantLIB version " << QL_VERSION << std::endl;
}

#endif /* check_version_hpp */
