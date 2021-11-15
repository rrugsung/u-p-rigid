/*************************************************************************
                        Material Point Method
                       Author: Shyamini Kularathna
                         University of Cambridge

NOTE:

FILE: PropertyParse.hpp
**************************************************************************/
#ifndef MPM_MISC_PROPERTYPARSE_H
#define MPM_MISC_PROPERTYPARSE_H
                                                                          
// c++ header files
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> 
#include <map>

// boost header files
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace mpm {
  namespace misc {
    void PARSE_INPUT_PARAMETERS(std::string& iLine);

    bool gravity;
    unsigned forceType;
    bool freeSurface;
    double dt_;
    double total_time_;
    double load_step_;
    unsigned time_factor_;
    unsigned scalar_beta_;
    double damping_coefficient_;
    double Gamma;
    double Beta;
    double permeability_;
    double compressibility_;
    double soil_depth_;
    std::map<std::string, double> propertyList;

    void READ_PROPERTIES(std::string& line);
  }
}

#include "PropertyParse.ipp"

#endif
