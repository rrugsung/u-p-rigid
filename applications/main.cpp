/*************************************************************************
                        2D-Material Point Method
                       Author: Shyamini Kularathna
                         University of Cambridge
Description: This file contains the incompressible fluid flow solver based on
             the Chorin's Projection method. The incompressible Navier-Stokes
             equations are solved using a semi-implicit approach in three steps.
             There are weveral options to be chosen by the user.
                  1. Gravity is known at nodes or particles
                  2. Consistent mass matrix or lumped mass matrix
**************************************************************************/
// c++ header files
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>

#include "file_handle.hpp"
#include "PropertyParse.hpp"
#include "Mesh.hpp"
#include "MpmParticle.hpp"
#include "Element.hpp"
#include "Node.hpp"
#include "Particle.hpp"
#include "projection_solver.hpp"                                                                        

int main (int argc, char* argv[]) {
  std::cout << "\n \n \t \t Semi-implicit Two Phase Material Point Method \n";
  std::cout << "\t \t \t Univerisy Of Cambridge \n";

  // check the existence of input folder
  if (argc != 2)
    std::cerr << "ERROR: in arguments" << "\n";

  boost::filesystem::path p (argv[1]);
  if (boost::filesystem::exists(p))
    {
      if (boost::filesystem::is_directory(p))
	std::cout << "\n \t Input data is from " << p << "\n";
      else
	std::cerr << "ERROR: " << p << " is not a directory \n";
    }
  else
    std::cerr << "ERROR: " << p << " does not exist" << "\n";

  mpm::FileHandle file_handle(p);
  mpm::Mesh* mesh = file_handle.read_mesh();
  mesh->iterate_over_elements(std::bind(&mpm::Element::compute_centre_coordinates, std::placeholders::_1));
  mpm::MpmParticle* particles = file_handle.read_particles();
  double const_permeability = mpm::misc::permeability_;
  double fluid_compressibility = mpm::misc::compressibility_;
  mpm::ProjectionSolver* solver = new mpm::ProjectionSolver(const_permeability, fluid_compressibility);

  double dt = mpm::misc::dt_;
  double total_analysis_time = mpm::misc::total_time_;
  double load_step = mpm::misc::load_step_;
  unsigned time_factor = mpm::misc::time_factor_;
  const double free_pressure = 0.; // prescribed pressure at the free surface
  const double material_density = 1000; // this should be improved
  std::cout << " Read All InputFiles" << "\n \n";
  
  unsigned write_steps = 0;
  long double total_time = 0.;
  long double sub_time = 0.;
  double analysis_time = total_analysis_time * time_factor;
  double accumulate_time = 0.;
  dt = dt * time_factor;
  double damping_factor = 0.;
  damping_factor = mpm::misc::damping_coefficient_;

  bool TwoPhase = false;
  if (particles->number_of_tp_particles() != 0){
    TwoPhase = true;
    std::cout << " Two-Phase Particles Included" << "\n \n";
  }

  bool Contact = true;

  //! MOVING MESH PARAMETERS
  double soil_depth = 4.0;
  double rigid_displacement;

  //! TIME STEP BEGINS
  unsigned step = 0;
  while (accumulate_time < total_analysis_time) {

    //! INITIALISE MESH AND PARTICLES
    mesh->initialise_mesh();
    particles->iterate_over_particles(std::bind(&mpm::Particle::initialise_particle, std::placeholders::_1));

    const auto begin = std::chrono::high_resolution_clock::now();
    //! WRITE OUTPUTS
    if (((load_step * time_factor * write_steps)- (accumulate_time * time_factor)) < dt * 0.000001) {
      std::cout << "\nStep: " << step << " " << accumulate_time
		<< "/" << total_analysis_time << " "
		<< sub_time << " ms "
		<< "total time: " << total_time << " ms \n";
      file_handle.write_data(write_steps, particles, mesh);
      write_steps++;
      sub_time = 0.;
      //step++;
    }  

    //! LOCATE PARTICLES IN THE MESH AND COMPUTE SHAPE FUNCTIONS
    mesh->locate_particles_in_mesh(particles);
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_local_coordinates, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_shape_functions, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_global_derivatives_shape_functions, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_global_derivatives_shape_functions_at_centre, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_B_matrix, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_BBar_matrix, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_B_matrix_at_centre, std::placeholders::_1));

    //! DEFINE FREE SURFACE
    if(mpm::misc::freeSurface) {
      mesh->iterate_over_elements_of_p(std::bind(&mpm::Element::find_free_surface_nodes, std::placeholders::_1));
      mesh->prescribe_pressure_at_free_surface(free_pressure);
    }

    //! FIND PARTIALLY FILLED ELEMENTS
    mesh->find_partially_filled_elements();
   
    //if (step == 1) {
      //if (TwoPhase) {
        //particles->iterate_over_tp_particles(std::bind(&mpm::Particle::map_pore_pressure_to_nodes, std::placeholders::_1));
        //mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_nodal_pressure, std::placeholders::_1));
        //particles->iterate_over_tp_particles(std::bind(&mpm::Particle::map_pore_pressure_from_nodes, std::placeholders::_1));
      //}
    //}

    //! MAP PARTICLE INFORMATION TO GRID NODES
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::update_nodal_phase, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::append_material_id_to_nodes, std::placeholders::_1));
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::map_mass_to_nodes, std::placeholders::_1, Contact));
    particles->iterate_over_sp_particles(std::bind(&mpm::Particle::map_sp_mass_to_nodes, std::placeholders::_1, Contact));
    particles->iterate_over_particles(std::bind(&mpm::Particle::map_momentum_to_nodes, std::placeholders::_1, Contact));
    mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_nodal_initial_velocity, std::placeholders::_1, Contact));
    
    //! COMPUTE STRAIN AND STRESSES
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_solid_strain_rate, std::placeholders::_1, Contact));
    particles->iterate_over_particles(std::bind(&mpm::Particle::compute_solid_strain, std::placeholders::_1,dt));
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::compute_solid_stress, std::placeholders::_1,dt));

    //! ASSIGN FORCES TO NODES
    //particles->iterate_over_tp_particles(std::bind(&mpm::Particle::assign_body_force_to_nodes, std::placeholders::_1, Contact));
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::assign_traction_force_to_nodes, std::placeholders::_1,Contact, accumulate_time));
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::assign_internal_force_to_nodes, std::placeholders::_1, Contact));
    
    particles->iterate_over_sp_particles(std::bind(&mpm::Particle::assign_sp_body_force_to_nodes, std::placeholders::_1, Contact));
    particles->iterate_over_sp_particles(std::bind(&mpm::Particle::assign_sp_traction_force_to_nodes, std::placeholders::_1, Contact, accumulate_time));
    particles->iterate_over_sp_particles(std::bind(&mpm::Particle::assign_sp_internal_force_to_nodes, std::placeholders::_1, Contact));

    // COMPUTE MULTIMATERIAL UNIT NORMAL VECTOR
    particles->iterate_over_particles(std::bind(&mpm::Particle::map_multimaterial_domain_gradients, std::placeholders::_1));
    mesh -> iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_multimaterial_normal_unit_vectors, std::placeholders::_1));
    particles->iterate_over_sp_particles(std::bind(&mpm::Particle::compute_penalty_factor, std::placeholders::_1));
    // mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_nodal_damping_forces, std::placeholders::_1,damping_factor));

    
    mesh->compute_rigid_body_int_acceleration();
    
    // COMPUTE INTERMEDIATE RIGID BODY ACCELERATION
    // BUILD MATRICES FOR THE CG SEMI-IMPLICIT SOLVER
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::compute_element_matrix_mp, std::placeholders::_1, dt));
    solver->assemble_solver(mesh);
    
    // COMPUTE INTERMEDIATE SOLID AND WATER VELOCITIES
    // i=0 for CoM and rigid body, i=1 for the first material (mat_id = 0).
    for (unsigned i=0; i<2; i++) {
      mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_intermediate_solid_acceleration, std::placeholders::_1, i, dt));
      mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_intermediate_water_velocity, std::placeholders::_1, i));
    }

    //! SOLVE PROJECTION METHOD: Compute acceleration, velocity, and pore pressure.
    solver->solve_pressure_poisson_equation(dt);
    solver->assign_final_pressure_to_nodes();
    solver->assign_final_acceleration_force_to_nodes();

    // COMPUTE FINAL ACCELERATIONS AND VELOCITIES
    mesh->compute_rigid_body_final_acceleration();
    mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_final_solid_acceleration, std::placeholders::_1, 1, dt));

    // APPLY CONTACT MECHANICS
    mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_multimaterial_relative_velocities, std::placeholders::_1));
    mesh -> iterate_over_nodes_of_p(std::bind(&mpm::Node::apply_contact_mechanics, std::placeholders::_1, dt));
    //mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::matrix_test, std::placeholders::_1));

    //! UPDATE PARTICLES
    if (Contact)
      particles->iterate_over_particles(std::bind(&mpm::Particle::update_contact_velocity_and_position, std::placeholders::_1, dt));
    else
      particles->iterate_over_particles(std::bind(&mpm::Particle::update_velocity_and_position, std::placeholders::_1, dt, accumulate_time));
        
    particles->iterate_over_tp_particles(std::bind(&mpm::Particle::update_pressure, std::placeholders::_1));
    particles->iterate_over_particles(std::bind(&mpm::Particle::update_porosity, std::placeholders::_1, dt));

    //! COMPUTE STRAIN AND STRESSES
    //particles->iterate_over_particles(std::bind(&mpm::Particle::compute_solid_strain_rate, std::placeholders::_1, Contact));
    //particles->iterate_over_particles(std::bind(&mpm::Particle::compute_solid_strain, std::placeholders::_1,dt));
    //particles->iterate_over_tp_particles(std::bind(&mpm::Particle::compute_solid_stress, std::placeholders::_1,dt));

    //! PRESSURE SMOOTHENING - DEACTIVATE IF NOT NECESSARY
    if (TwoPhase){
      particles->iterate_over_tp_particles(std::bind(&mpm::Particle::map_pore_pressure_to_nodes, std::placeholders::_1));
      mesh->iterate_over_nodes_of_p(std::bind(&mpm::Node::compute_nodal_pressure, std::placeholders::_1));
      particles->iterate_over_tp_particles(std::bind(&mpm::Particle::map_pore_pressure_from_nodes, std::placeholders::_1));
    }

    rigid_displacement = particles->give_rigid_displ();
    mesh->iterate_over_nodes(std::bind(&mpm::Node::update_mesh_configuration, std::placeholders::_1, rigid_displacement, soil_depth));
    soil_depth = soil_depth + rigid_displacement;

    auto step_time = std::chrono::high_resolution_clock::now() - begin;
    auto duration = std::chrono::duration <double, std::milli> (step_time).count();
    total_time += duration;
    sub_time += duration;
    accumulate_time += (dt / time_factor);
    step++;
  }
   delete mesh;
   delete particles;
   delete solver;
   return 0;
}

