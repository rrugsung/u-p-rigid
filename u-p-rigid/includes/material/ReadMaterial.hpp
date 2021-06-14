/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: ReadMaterial.hpp
****************************************************************************/

#ifndef MPM_MATERIAL_READMATERIAL_H
#define MPM_MATERIAL_READMATERIAL_H

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// c++ header files
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

// mpm header files
#include "MaterialBase.hpp"
#include "MohrCoulomb/MohrCoulomb.hpp"
#include "ILE/ILE.hpp"
#include "PropertyParse.hpp"

namespace mpm {
    namespace material {
        class ReadMaterial;
    }
}

class mpm::material::ReadMaterial {

protected:
    typedef mpm::material::MaterialBase* MaterialBasePtr;
    typedef std::vector<MaterialBasePtr> MaterialVector;


public:
  // CONSTRUCTOR
  //! param[in] matFile input file containing material model parameters
  ReadMaterial(std::ifstream& matFile);

private:
  // register materials
  //! param[in] name material model name
  mpm::material::MaterialBase* registerMaterial(std::string& name);

public:
  // give pointers to the registered material models
  //! param[out] materialPtrs_ vector of pointers to the registered materials
  std::vector<mpm::material::MaterialBase*> givePtrsToMaterials() {
    return materialPtrs_;
  }
  
  unsigned numMatTypes;
  
public:
  std::vector<mpm::material::MaterialBase*> materialPtrs_;
};

#include "ReadMaterial.ipp"

#endif
