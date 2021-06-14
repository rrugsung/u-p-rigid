
#ifndef MPM_MATERIAL_MATERIALBASE_H
#define MPM_MATERIAL_MATERIALBASE_H

#include <iostream>
#include <vector>
#include <map>

#include "PropertyParse.hpp"
#include "Constants.hpp"

namespace mpm {
    namespace material {
        class MaterialBase;
    }
}

class mpm::material::MaterialBase {

protected:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned dof = 3*(dim-1);

public:
  virtual double give_density() = 0;
  virtual double give_porosity() = 0;
  virtual double give_permeability() = 0;
  virtual void compute_stress(const Eigen::Matrix<double,1,dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& pressure) = 0;
  virtual void change_friction_coefficient() = 0;

  void setProperty(std::string propName, double& propValue) {
    propValue = mpm::misc::propertyList[propName];
  }

};

#endif
