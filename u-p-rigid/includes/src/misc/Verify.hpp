/*************************************************************************
                        Material Point Method
                       Author: Shyamini Kularathna
                         University of Cambridge
NOTE: This code is designed 

FILE: Verify.hpp
**************************************************************************/
#ifndef MPM_MISC_VERIFY_H
#define MPM_MISC_VERIFY_H

// c++ header files
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream> 

namespace mpm {
    namespace misc {
        void VERIFY_OPEN(std::ifstream& file, std::string& fileName);
    }

}

#include "Verify.ipp"

#endif
