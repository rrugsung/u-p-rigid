/*************************************************************************
                        Material Point Method
                       Author: Shyamini Kularathna
                         University of Cambridge
FILE: Projection.hpp
**************************************************************************/
#ifndef MPM_PROJECTION_SOLVER_H
#define MPM_PROJECTION_SOLVER_H 
                                                                         
// c++ header files
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <tuple>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "Constants.hpp"
#include "Mesh.hpp"
#include "Node.hpp"
#include "Element.hpp"
#include "element_matrix.hpp"
#include "MpmParticle.hpp"


namespace mpm {
    class ProjectionSolver;
  template <class T>
  struct setcomparator {
    bool operator () (const T& arg1, const T& arg2) const {
      unsigned firstId = arg1 -> give_id();
      unsigned secondId = arg2 -> give_id();
      return firstId < secondId;
    }
  };
}


class mpm::ProjectionSolver {
protected:
    static const unsigned dim = mpm::constants::DIM;
    static const unsigned numNodes = mpm::constants::NUMNODES;
    static constexpr double theta = 0.28;
  static constexpr double permeability_ = 0.0;

public:
  // constructor
  ProjectionSolver(const double &permeability, const double &compressibility);

  // initialise the system
  //! param[in] mesh_ptr pointer to the mesh
  void assemble_solver(mpm::Mesh* mesh_ptr);

  // solve pressure poisson equation
  //! param[in] dt time step
  //! param[in] mat_density material density
  bool solve_pressure_poisson_equation(const double &dt);

  // solve final velocity
  //! param[in] dt time step
  //! param[in] mat_density material density 
  bool solve_final_acceleration(const double &dt, const double & gamma);

  // assign solved pressure coefficients to nodes
  void assign_final_pressure_to_nodes();

  // assign solved final divergence free velocity to nodes
  void assign_final_acceleration_force_to_nodes();

  // give number of nodes for CG solver
  unsigned give_CG_nodes_size() const {
    return num_nodes_;
  }

  

  void print_nnodes();

private:
    // evaluate element mapping matrix
    void element_mapping_matrix();

    // store pressure constraints at nodes
    void store_pressure_constraints();

    // apply pressure boundary conditions to the linear system
    void apply_pressure_boundary_conditions_to_system();
   
    // apply velocity boundary conditions to the linear system
    void apply_velocity_boundary_conditions_to_system();

    // apply solid velocity boundary conditions to the linear system
    void apply_solid_velocity_boundary_conditions_to_system();

    // apply water velocity boundary conditions to the linear system
    void apply_water_velocity_boundary_conditions_to_system();

    // build global mass matrix
    void build_M_matrix();

    // build global K bar matrix
    void build_Kbar_matrix();

    // build global laplace matrix
    void build_L_matrix();

    // build global K_L2 and K_L3 matrices
    void build_K_L_matrix();

    // build global K_L4 matrix
    void build_K_L4_matrix();

    // build K_S2 matrix
    void build_K_S2_matrix();

    // build K_W2 matrix
    void build_K_W2_matrix();

    // check for zeroes in the stiffnes matrix and force vector
    void check_for_zeroes_in_linear_system();

    // solve linear equation using Conjugate Gradient method
    //! param[in] force_vector RHS force vector
    //! param[in] stiffnes_matrix LHS stiffness matrix
    //! param[in/out] solution 
    bool Conjugate_Gradient(Eigen::VectorXd &solution);
    bool Least_Squares_Conjugate_Gradient(Eigen::VectorXd &solution);

protected:
  // number of total nodes in the domain
  unsigned num_nodes_;
  // velocity dogree of freedom
  unsigned velocity_dof_;
  // pressure degree of freedom
  unsigned pressure_dof_;

  // vector of pointers to elements in the material domain
  std::vector<mpm::Element*> element_ptrs_;
  // vector of pointers to nodes in the material domain
  std::vector<mpm::Node*> node_ptrs_;

  // pressure constraints
  std::vector<std::tuple<bool,double>> pressure_constraints_;
  // velocity constraints
  std::vector<std::tuple<bool,double>> velocity_constraints_;
  // solid velocity constraints
  std::vector<std::tuple<bool,double>> solid_velocity_constraints_;
  // water velocity constraints
  std::vector<std::tuple<bool,double>> water_velocity_constraints_;

  // element mapping matrix
  Eigen::MatrixXi element_map_;

  // global mass matrix
  Eigen::SparseMatrix<double> M_;
  Eigen::SparseMatrix<double> M_11_;
  Eigen::SparseMatrix<double> M_22_;
  // global K Bar matrix
  Eigen::SparseMatrix<double> K_bar_;
  Eigen::SparseMatrix<double> K_single_;
  // global laplace matrix
  Eigen::SparseMatrix<double> L_;
  // global K_L2 matrix
  Eigen::SparseMatrix<double> K_L2_;
  // global K_L3 matrix
  Eigen::SparseMatrix<double> K_L3_;
  // global K_L4 matrix
  Eigen::SparseMatrix<double> K_L4_;
  // global K_S2 matrix
  Eigen::SparseMatrix<double> K_S2_;
  // global K_W2 matrix
  Eigen::SparseMatrix<double> K_W2_;


  // stiffness matrix of the linear equation
  Eigen::SparseMatrix<double> stiffness_matrix_;
  // force vector of the linear equation
  Eigen::VectorXd force_vector_;

  // initial velocity at nodes
  Eigen::VectorXd known_solid_velocity_;
  Eigen::VectorXd known_water_velocity_;
  // intermediate velocity at nodes
  Eigen::VectorXd intermediate_acceleration_;
  
  Eigen::VectorXd intermediate_solid_velocity_;
  Eigen::VectorXd intermediate_water_velocity_;
  // final velocity
  Eigen::VectorXd final_solid_acceleration_;
  Eigen::VectorXd final_water_acceleration_;
  // pressure projection coefficient
  Eigen::VectorXd pressure_;

  double const_permeability_;
  double fluid_compressibility_;
};

#include "projection_solver.ipp"

#endif

// dim      : 1D/2D/3D
// numNodes : number of nodes of the element
// theta    : cutoff ratio (element particle volume/element volume)
              //! elements below this cutoff value will not consider
              //! to be in material domain
