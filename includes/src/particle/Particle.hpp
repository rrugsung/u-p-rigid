/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         University Of Cambridge
File: Particle.hpp
****************************************************************************/
#ifndef MPM_PARTICLE_H
#define MPM_PARTICLE_H

#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>

#include <eigen3/Eigen/Dense>

#include "Constants.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "MaterialBase.hpp"

namespace mpm {
  class Particle;
}

class mpm::Particle {
protected:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static const unsigned dof = 3*(dim-1);
  unsigned beta_scalar_ = mpm::misc::scalar_beta_;
  static const unsigned numMats = mpm::constants::NUMBODIES;

public: 
  // constructor
  //! param[in] id particle id
  //! param[in] mat_id material id
  //! param[in] spacing particle initial spacing
  //! param[in] number of phase
  Particle(const unsigned& id, const unsigned& mat_id, const Eigen::Matrix<double,1,dim>& spacing,const double& vol, const unsigned& nphase);

  // initialise the particle
  void initialise_particle();

  // set element and nodes
  //! param[in] ptr_e pointer to the element
  //! param[in] ptrs_n vector of pointers to nodes of the element
  void set_element_and_nodes(mpm::Element* &ptr_e, Eigen::Matrix<mpm::Node* ,1, numNodes> &ptrs_n) {
        element_ = ptr_e;
        nodes_   = ptrs_n;
  }

  // set particle coordinates
  //! param[in] coordinates
  void set_coordinates(const Eigen::Matrix<double,1,dim> &coordinates) {
      coord_ = coordinates;
  }

  // print out particle id
  void print_id();

  // set particle initial pore pressure
  //! param[in] init_pressure initial pore pressure at particle
  void set_initial_pore_pressure(double &initial_pressure);

  // set particle initial stress
  //! param[in] init_stress initial stress vector at particle
  void set_initial_stress(const Eigen::Matrix<double, 1, 6> &initial_stress);

  // set solid particle surface traction
  //! param[in] direction direction of the surface
  //! param[in] solidtraction value of the solid surface traction
  void set_solid_surface_traction(const unsigned &direction, const double &solidtraction);

  // set water particle surface traction
  //! param[in] direction direction of the surface
  //! param[in] watertraction value of the water surface traction
  void set_water_surface_traction(const unsigned &direction, const double &watertraction);

  // set particle material
  //! param[in] materials vector of pointers to the all materials
  void set_material(const std::vector<mpm::material::MaterialBase*> materials);
  void update_material(const unsigned &id) {
    material_ = material_ptrs_.at(id);
  }

  double give_rigid_displacement() {
    return displacement_(1);
  }

//####################################
  // upadate nodal phase
  void update_nodal_phase();

  //! Assign material id of this particle to nodes
  void append_material_id_to_nodes();
  
  // map two-phase particle mass to nodes
  void map_mass_to_nodes(bool contact);

  // map multimaterial particle masses to nodes
  void map_multimaterial_masses_to_nodes();

  // map single-phase particle mass to nodes
  void map_sp_mass_to_nodes(bool contact);

  // map particle volume to nodes : delete if not in use
  void map_volume_to_nodes();

  // map particle momentum to nodes
  void map_momentum_to_nodes(bool contact);

  // map particle acceleration to nodes
  void map_acceleration_to_nodes();

  // map particle porosity to nodes
  void map_porosity_to_nodes();

  // map particle pressure to nodes
  void map_pore_pressure_to_nodes();

  // map pressure from nodes
  void map_pore_pressure_from_nodes();
  void map_solid_velocity_from_nodes();

//#####################################

  // assign two-phase body force to nodes
  void assign_body_force_to_nodes(bool contact);

  // assign single-phase body force to nodes
  void assign_sp_body_force_to_nodes(bool contact);

  // assign two-phase traction force to nodes
  void assign_traction_force_to_nodes(bool contact, const double &time);
  
  // assign single-phase traction force to nodes
  void assign_sp_traction_force_to_nodes(bool contact, const double &time);

  // assign two-phase internal force to nodes
  void assign_internal_force_to_nodes(bool contact);

  // assign single-phase internal force to nodes
  void assign_sp_internal_force_to_nodes(bool contact);

//#####################################
  void map_multimaterial_domain_gradients();

  void compute_penalty_factor();
//#####################################

  // compute strain rate
  void compute_solid_strain_rate(bool contact);

  // compute strain rate at element center
  void compute_solid_strain_rate_at_elem_center(bool contact);

  // compute strain
  void compute_solid_strain(const double dt);

  // compute stress
  void compute_solid_stress(const double dt);

//#####################################
  // update interface velocity
  void update_contact_velocity_and_position(const double& dt);
  
  // update velocity
  void update_velocity_and_position(const double& dt, const double& time);

  // update position
  void update_position(const double& dt);

  // update pressure
  void update_pressure();

  // update density
  void update_density_and_mass(const double &compressibility);
 
  // update porosity
  void update_porosity(const double &dt);

//####################################

  // compute local coordinates
  void compute_local_coordinates();

  // compute shape functions
  void compute_shape_functions();

  // compute global derives of shape functions
  void compute_global_derivatives_shape_functions();

  // compute global derivatives of shape function at element centre
  void compute_global_derivatives_shape_functions_at_centre();

  // compute B matrix
  //!           B_ = [B_0, B_1, B_2, B_3]
  //!           where B_i = [B_ix    0      0 ]
  //!                       [ 0     B_iy    0 ]
  //!                       [ 0      0    B_iz]
  //!                        ------------------
  //!                       [B_iy   B_ix    0 ]
  //!                       [ 0     B_iz  B_iy]
  //!                       [B_iz    0    B_ix]
  //!            where B_ij = dN_i(x_p) / dj; j = x,y,z           
  //!
  void compute_B_matrix();
  void compute_BBar_matrix();

  // compute B matrx at the centre of the element
  void compute_B_matrix_at_centre();

  // build element matrix
  void compute_element_matrix_mp(const double &dt);

//###########################################

  // give particle id
  //! param[out] id_ particle id
  unsigned give_id() const {
    return id_;
  }

  // give particle coordinates
  //! param[out] coord_ particle coordinates
  Eigen::Matrix<double, 1, dim> give_coordinates() const {
    return coord_;
  }

  // give particle spacing
  //! param[out] spacing_ particle initial spacing
  Eigen::Matrix<double, 1, dim> give_spacing() const {
    return spacing_;
  }

  // give particle volume
  //! param[out] volume_ particle volume
  double give_volume() const {
    return volume_;
  }

  double give_density() const {
    return solid_grain_density_;
  }

  // give particle material id
  //! param[out] mat_id particle material id
  unsigned give_mat_id() const {
    return mat_id_;
  }

  // give particle phase
  //! param[out] nphase number of particle phase
  unsigned give_nphase() const {
    return nphase_;
  }
//###############################################

  // WRITE VELOCITY TO FILE
  void write_velocity(std::ostream& oFile);

  // WRITE PRESSURE TO FILE
  void write_pressure(std::ostream& oFile);

  // WRITE STRESS TO FILE
  void write_stress(std::ostream& oFile);

  // WRITE STRAIN TO FILE
  void write_strain(std::ostream& oFile);

  // WRITE DISPLACEMENT TO FILE
  void write_displacement(std::ostream& oFile);

  // TEMPORARY: WRITE PRESSURE TO DAT FILE
  void write_deviatoric_shear_strain(std::ofstream& oFile);
  void write_equivalent_plastic_strain(std::ofstream& oFile);

  // ###############################################
  // TEMPORARY
  void assign_moving_undrained_boundary() {
    if(id_ > 17879 && id_ < 17910) {
      for(unsigned i = 0; i < numNodes; i++)
	nodes_(i)->set_moving_undrained_boundary();
    }
  }
  // ###############################################

private:

  // compute particle mass
  void compute_mass();

  // change the sign of value
  void sign(double& variable, double value);

protected:
  unsigned id_;
  Eigen::Matrix<mpm::Node*,1, numNodes> nodes_;
  mpm::Element* element_;

  Eigen::Matrix<double, 1, dim> coord_;
  Eigen::Matrix<double, 1, dim> spacing_;
  double volume_;
  Eigen::Matrix<double, 1, dim> displacement_;
  unsigned nphase_;

  unsigned mat_id_;
  std::vector<mpm::material::MaterialBase*> material_ptrs_;
  mpm::material::MaterialBase* material_;
  Eigen::Matrix<double, 1, dim> gravity_;

  double porosity_;
  double permeability_;

  // surface traction
  std::vector<std::tuple<unsigned, double> > solid_traction_;
  std::vector<std::tuple<unsigned, double> > water_traction_;

  double solid_grain_density_;
  double water_grain_density_;
  double solid_mass_;
  double water_mass_;

  // Acceleration and Velocity
  Eigen::Matrix<double, 1, dim> solid_velocity_;
  Eigen::Matrix<double, 1, dim> water_velocity_;
  Eigen::Matrix<double, 1, dim> solid_acceleration_;
  Eigen::Matrix<double, 1, dim> water_acceleration_;

  double pore_pressure_;
  Eigen::Matrix<double, 1, 6> solid_stress_; // effective stress
  Eigen::Matrix<double, 1, dof> solid_strain_;
  Eigen::Matrix<double, 1, dof> deviatoric_strain_;
  Eigen::Matrix<double, 1, dof> solid_strain_rate_;
  Eigen::Matrix<double, 1, dof> solid_centre_strain_rate_;

  double plastic_strain_;
  double dev_shear_strain_;
  // double vol_strain_rate_;
  // double centre_vol_strain_rate_;

  Eigen::Matrix<double, 1, dim> xi_; // local coordinates
  Eigen::Matrix<double, 1, numNodes> shape_fun_;
  Eigen::Matrix<double, 1, numNodes> shape_fun_centre_;
  Eigen::Matrix<double, dim, numNodes> grad_shape_fun_;
  Eigen::Matrix<double, dim, numNodes> grad_shape_fun_centre_;

  std::array<Eigen::Matrix<double, dof, dim>, numNodes> B_;
  std::array<Eigen::Matrix<double, dof, dim>, numNodes> BBar_;
  std::array<Eigen::Matrix<double, dof, dim>, numNodes> BCentre_;

  Eigen::Matrix<double,1,dof> m;
};

#include "Particle.ipp"

#endif
