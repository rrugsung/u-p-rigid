/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: MohrCoulomb.hpp
****************************************************************************/

#ifndef MPM_MATERIAL_MOHRCOULOMB_H
#define MPM_MATERIAL_MOHRCOULOMB_H

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// c++ header files
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/math/constants/constants.hpp>

// mpm header files
#include "MaterialBase.hpp"
#include "PropertyParse.hpp"
#include "Constants.hpp" 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

namespace mpm {
    namespace material {
        class MohrCoulomb;
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class mpm::material::MohrCoulomb : public mpm::material::MaterialBase {
  const double pi = boost::math::constants::pi<double>();
  const double onethirdpi = 1.047197551;
public:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static const unsigned dof = 3 * (dim - 1);

  typedef Eigen::Matrix<double, dof,1> VectorDDOF;
  typedef Eigen::Matrix<double, 6,1>   VectorD6x1;
  typedef Eigen::Matrix<double, 1,6>   VectorD1x6;
  typedef Eigen::Matrix<double, 6, 6> Matrix6x6;

public:
  // CONSTRUCTOR
  MohrCoulomb();

  //! CREATE THE MATERIAL
  static MaterialBase* create() {
    return new MohrCoulomb();
  }

  //! GIVE DENSITY
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

  //! compute elastic stiffness matrix
  void compute_elastic_stiffness_matrix();

  void change_friction_coefficient() {}

  //! COMPUTE STRESS
  //! Mohr Coulomb model is implemented based on work by Dilan ROBERT
  // param[in] dstrain: strain increment
  // param[in] stress: stress
  // param[in] epds: equivalent plastic deviatoric strain
  // param[out] stress: stress is updated
  // param[out] epds: equivalent plastic deviatoric strain is updated
  void compute_stress(const Eigen::Matrix<double,1,dof>& dstrain,
		      Eigen::Matrix<double,1,6>& stress, double& epds);

private:  
  void compute_rho_theta(const VectorD6x1 stress, double& _j2, double& _j3);
  Eigen::Matrix<double,6,2> compute_dFdP(const double _j2,const double _j3,
					 const VectorD6x1 _devstress,
					 double& _psoftening);
  void compute_lambda_trial();


protected:
  Matrix6x6 De_;
  double density_;
  double porosity_;
  double permeability_;
  double E_;
  double mu_;
  double t_;

  // initial values
  double phi_;
  double c_;
  double psi_;

  // residual values
  double phi_resd_;
  double c_resd_;
  double psi_resd_;

  // current values
  double phi_cur_;
  double c_cur_;
  double psi_cur_;

  // plastic deviatoric strain
  double epds_peak_;
  double epds_crit_;

  // plastic deviatori strain vector
  double epds_;
  Eigen::Matrix<double,6,1> PDS_;

  // parameters for yield function
  double rho_;
  double theta_;
  double epsilon_;
};

#include "MohrCoulomb.ipp"

#endif
