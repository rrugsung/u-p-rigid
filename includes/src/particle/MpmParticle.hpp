/*****************************************************************************
                        Material Point Method
                         Shyamini Kularathna
                       University of Cambridge
FILE: MpmParticle.hpp
*****************************************************************************/
#ifndef MPM_MPMPARTICLE_H
#define MPM_MPMPARTICLE_H

// c++ header files
#include <vector>
#include <fstream>
#include <iostream>

#include <eigen3/Eigen/Dense>

#include "Constants.hpp"
#include "Particle.hpp"
#include "MaterialBase.hpp"

namespace mpm {
  class MpmParticle;
}

class mpm::MpmParticle {

protected:
    static const unsigned dim = mpm::constants::DIM;
    static const unsigned numNodes = mpm::constants::NUMNODES;
    static const unsigned dof = 3*(dim-1);

public:
  // constructor
  MpmParticle();

  // destructor
  //! do nothing
  ~MpmParticle() { }

  // read particles
  //! param[in] p_file input file for particle data
  //! param[in] s_file input file for initial stress
  void read_particles(std::ifstream& p_file, std::ifstream& s_file);

  // Read traction at particles
  //! param[in] tractionFile input file for traction force at particles
  void read_surface_traction(std::ifstream& traction_file);

  // assign material to particle
  //! param[in] material_ptrs vector of pointers to all materials
  void assign_material_to_particles(std::vector<mpm::material::MaterialBase*> &material_ptrs);

  // Iterate over particles
  template<typename FP>
  void iterate_over_particles(FP function);

  // Iterate over particles
  template<typename FP>
  void iterate_over_sp_particles(FP function);

  // Iterate over particles
  template<typename FP>
  void iterate_over_tp_particles(FP function);

  // give total number of particles
  // ! param[out] number of particles
  unsigned number_of_particles() {
    return particles_.size();
  }

  // give number of two-phase particles
  unsigned number_of_tp_particles(){
    return tp_particles_.size();
  }
  
  // give pointer to the particle of given id
  //! param[out] pointer to the particle at given id
  std::shared_ptr<mpm::Particle> pointer_to_particle(const unsigned& id) {
    return particles_.at(id);
  }

  double give_rigid_displ ();
  
  // give coordinates of the particle of given id
  //! param[out] coordinates of the particle at given id
  Eigen::Matrix<double, 1, dim> particle_coordinates(const unsigned& id) {
    Eigen::Matrix<double,1,dim> p_coord = particles_.at(id)->give_coordinates();
    return p_coord;
  }

  // WRITE PARTICLE VELOCITY DATA
  void write_particle_velocity_data_to_file(std::ostream& outFile);

  // WRITE PARTICLE PRESSURE DATA
  void write_particle_pressure_data_to_file(std::ostream& outFile);

  // WRITE PARTICLE STRESS DATA
  void write_particle_stress_data_to_file(std::ostream& outFile);

  // WRITE PARTICLE STRAIN DATA
  void write_particle_strain_data_to_file(std::ostream& outFile);

  // WRITE PARTICLE DISPLACEMENT DATA
  void write_particle_displacement_data_to_file(std::ostream& outFile);

  // temporary function
  void write_particle_deviatoric_strain_data_to_file(std::ofstream &outFile);
  void write_particle_plastic_strain_data_to_file(std::ofstream &outFile);

protected:
  // Pointers to the all particles
  std::vector<std::shared_ptr<mpm::Particle>> particles_;

    // Pointers to the two-phase particles
  std::vector<std::shared_ptr<mpm::Particle>> tp_particles_;

    // Pointers to the single-phase particles
  std::vector<std::shared_ptr<mpm::Particle>> sp_particles_;
};

#include "MpmParticle.ipp"

#endif
