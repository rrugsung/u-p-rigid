/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: MohrCoulomb.hpp
****************************************************************************/

#ifndef MPM_MATERIAL_MOHRCOULOMB_H
#define MPM_MATERIAL_MOHRCOULOMB_H

// c++ header files
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/math/constants/constants.hpp>

// mpm header files
#include "MaterialBase.hpp"
#include "Constants.hpp"

namespace mpm {
    namespace material {
        class MohrCoulomb;
    }
}

class mpm::material::MohrCoulomb : public mpm::material::MaterialBase {
  const double pi = boost::math::constants::pi<double>();
protected:
    static const unsigned dim = mpm::constants::DIM;
    static const unsigned dof = 3*(dim-1);

public:
  // CONSTRUCTOR
  MohrCoulomb();

  // create the material
  //! param[out] pointer to the material
  static MaterialBase* create() {
    return new MohrCoulomb();
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

  // COMPUTE STRESS
  void compute_stress(const Eigen::Matrix<double, 1, dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& plastic_strain);

  void compute_elastic_stiffness_matrix();
  void change_friction_coefficient() {
    phi_ = 25;
    //c_= 50000.;
  }

protected:
  double density_;
  double porosity_;
  double permeability_;
  double E_;
  double mu_;
  double phi_;
  double c_;
  double psi_;
  double t_;

  Eigen::Matrix<double,6,6> De_;

  double phi_crit_;
  double c_crit_;
  double gamma_peak_;
  double gamma_crit_;

  Eigen::Matrix<double, 1, dof> pstrain_;
};

#include "MohrCoulomb.ipp"

#endif
