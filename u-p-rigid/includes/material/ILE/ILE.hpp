/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: ILE.hpp
****************************************************************************/

#ifndef MPM_MATERIAL_ILE_H
#define MPM_MATERIAL_ILE_H

// c++ header files
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

// mpm header files
#include "MaterialBase.hpp"
#include "PropertyParse.hpp"
#include "Constants.hpp" 

namespace mpm {
    namespace material {
        class ILE;
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class mpm::material::ILE : public mpm::material::MaterialBase {

protected:
    static const unsigned dim = mpm::constants::DIM;
    static const unsigned numNodes = mpm::constants::NUMNODES;
    static const unsigned dof = 3 * (dim - 1);

public:
  // CONSTRUCTOR
  ILE();

  // create the material
  //! param[out] pointer to the material
  static MaterialBase* create() {
    return new ILE();
  }

  // give density
  //! param[out] density_
  double give_density() {
    return density_;
  }

  // give porosity
  //! param[out] porosity_
  double give_porosity() {
    return porosity_;
  }

  // give permeability
  //! param[out] permeability_
  double give_permeability() {
    return permeability_;
  }

  // compute stress
  void compute_stress(const Eigen::Matrix<double, 1, dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& pressure);

  void compute_elastic_stiffness_matrix();
  void change_friction_coefficient() {}

protected:
  double density_;
  double porosity_;
  double permeability_;
  double E_;
  double mu_;
  Eigen::Matrix<double, 6, 6> De_;
};

#include "ILE.ipp"

#endif
