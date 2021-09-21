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
#include "Material_utility.hpp"

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

  enum FailureState {Elastic, Tensile, Shear};

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

  //! Compute yield function and yield state
  //! \param[in] state_vars History-dependent state variables
  //! \retval yield_type Yield type (elastic, shear or tensile)
  mpm::material::MohrCoulomb::FailureState compute_yield_state(
    Eigen::Matrix<double, 2, 1>* yield_function, double epsilon, double rho, double theta, double phi, double cohesion);

private:  
  void compute_stress_invariants(const VectorD6x1 stress);
  void compute_df_dp(mpm::material::MohrCoulomb::FailureState yield_type,
    const VectorD6x1& stress, VectorD6x1* df_dsigma, VectorD6x1* dp_dsigma, const double rho, const double theta, 
    const double phi, const double psi, double pdstrain, double* dp_dq, double* softening);
  void compute_lambda_trial();


protected:
  Matrix6x6 de_;
  double density_;
  double porosity_;
  double permeability_;

  double phi_cur;
  double psi_cur;
  double cohesion_cur;

  // plastic deviatori strain vector
  double epds_;
  Eigen::Matrix<double,6,1> PDS_;

  // parameters for yield function
  double rho_;
  double theta_;
  double epsilon_;

  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Shear modulus
  double shear_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Maximum friction angle phi
  double phi_peak_{std::numeric_limits<double>::max()};
  //! Maximum dilation angle psi
  double psi_peak_{std::numeric_limits<double>::max()};
  //! Maximum cohesion
  double cohesion_peak_{std::numeric_limits<double>::max()};
  //! Residual friction angle phi
  double phi_residual_{std::numeric_limits<double>::max()};
  //! Residual dilation angle psi
  double psi_residual_{std::numeric_limits<double>::max()};
  //! Residual cohesion
  double cohesion_residual_{std::numeric_limits<double>::max()};
  //! Peak plastic deviatoric strain
  double pdstrain_peak_{std::numeric_limits<double>::max()};
  //! Residual plastic deviatoric strain
  double pdstrain_residual_{std::numeric_limits<double>::max()};
  //! Tension cutoff
  double tension_cutoff_{std::numeric_limits<double>::max()};
};

#include "MohrCoulomb.ipp"

#endif
