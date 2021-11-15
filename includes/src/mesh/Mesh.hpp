/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         Unversity Of Cambridge
File: Mesh.hpp
****************************************************************************/

#ifndef MPM_MESH_H
#define MPM_MESH_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <set>
#include <functional>

#include <eigen3/Eigen/Dense>

#include "Constants.hpp"
#include "Element.hpp"
#include "MpmParticle.hpp"
#include "Node.hpp"

namespace mpm {
  class Mesh;
  template <class T>
  struct comparator {
    bool operator () (const T& arg1, const T& arg2) const {
      unsigned firstId = arg1 -> give_id();
      unsigned secondId = arg2 -> give_id();
      return firstId < secondId;
    }
  };
}

class mpm::Mesh {

protected:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static const unsigned numCorners = ((dim - 1) * dim) + 2;
  static constexpr double beta = 0.7;

public:
  // default constructor
  //! param[in] meshFile meshData input file
  Mesh(std::ifstream& mesh_data_file);

  // default destructor
  //! free all dynamically allocated memory
  ~Mesh() { free_memory(); }

  // read nodes and elements
  //! param[in] node_file input file for nodes
  //! param[in] elem_file input file for elements
  void read_nodes_and_elements(std::ifstream& node_file, std::ifstream& elem_file);

  // read general constraints
  //! param[in] gen_con input file for general constraints
  void read_general_constraints(std::ifstream& vel_con_file);

  // read frinction constraints
  //! param[in] fric_con input file for frictional constraints
  void read_friction_constraints(std::ifstream& fric_con_file);

  // locate particles in mesh
  //! param[in] particle_set vector of pointers to all particles
  void locate_particles_in_mesh(mpm::MpmParticle* &particle_set);
 
  // initialise mesh
  void initialise_mesh();

  // find partially filled elements
  //! param[in] ratio limiting ratio of particle volume/element volume
  void find_partially_filled_elements();

  // prescribe the pressure at nodes located at free surface
  //! param[in] pressure prescribed pressure value
  void prescribe_pressure_at_free_surface(const double &pressure);

  void compute_rigid_body_initial_velocity();

  void compute_rigid_body_int_acceleration(const double &time);

  void compute_rigid_body_final_acceleration(const double &time);

  template<typename FP>
  void iterate_over_elements(FP function) const;

  template<typename FP>
  void iterate_over_nodes(FP function) const;

  template<typename FP>
  void iterate_over_elements_of_p(FP function) const;

  template<typename FP>
  void iterate_over_nodes_of_p(FP function) const;

  // give number of elements which contain particles
  unsigned give_num_elements_of_p() {
    return p_element_set_.size();
  }

  // give number of nodes which contain particles
  unsigned give_num_nodes_of_p() {
    return p_node_set_.size();
  }

  // give vector of pointers to elements which contain particles
  //! param[in/out] vec_of_elem vector of pointers to elements
  void give_elements_of_p(std::vector<mpm::Element*>& vec_of_elem) {
    for (const auto &i : p_element_set_)
      vec_of_elem.push_back(i);
  }

  // give vector of pointers to elements which contain particles above cutoff
  //! param[in/out] vec_of_elem vector of pointers to elements
    void give_elements_with_particles_above_cutoff(std::vector<mpm::Element*>& vec_of_elem, const double &cut_off);

  // give vector of pointers to nodes which contain particles
  //! param[in/out] vec_of_nodes vector of pointers to nodes
  void give_nodes_of_p(std::vector<mpm::Node*>& vec_of_nodes) {
    for (const auto &i : p_node_set_)
      vec_of_nodes.push_back(i);
  }

  // WRITE MESH DATA
  void write_mesh_data_to_file(std::ostream& outFile);

  // WRITE PARA DATA
  void write_para_data_to_file(std::ostream& outFile);

private:
  // free all dynamically allocated memory within "MeshBase" class
  void free_memory();

  void check_particle_is_inside_mesh(Eigen::Matrix<int, 1 , dim> &eGrid, unsigned &pId);

  void set_elements_and_nodes_of_particles(unsigned &element_id, std::shared_ptr<mpm::Particle> &particle_ptr);

protected:
  // Pointers to all elements in the mesh
  std::vector<mpm::Element*> elements_;
  // Pointers to all nodes in the mesh
  std::vector<mpm::Node*> nodes_;

  // Set of elements which containes particles
  std::set<mpm::Element*, mpm::comparator<mpm::Element*> > p_element_set_;
  std::set<mpm::Element*, mpm::comparator<mpm::Element*> > solver_element_set_;
  // Set of nodes which containes particles
  std::set<mpm::Node*, mpm::comparator<mpm::Node*> > p_node_set_;
  // Set of partially filled elements
  std::vector<mpm::Element*> part_fill_elems_;
  // Set of nodes located at free surface
  std::set<mpm::Node*, mpm::comparator<mpm::Node*> > free_node_set_;


  // Mesh spacing
  Eigen::Matrix<double, 1, dim> mesh_spacing_;
  // number of oelements in each directions
  Eigen::Matrix<unsigned, 1, dim> num_elements_;
  // Ids of the corner elements of the mesh
  Eigen::Matrix<unsigned, 1, numCorners> corner_elements_;
  // Ids of the corner noes of the mesh
  Eigen::Matrix<unsigned, 1, numCorners> corner_nodes_;
  // coordinates of the first node of the mesh
  Eigen::Matrix<double, 1, dim> first_node_coord_;
  // coordinates of the last node of the mesh
  Eigen::Matrix<double, 1, dim> last_node_coord_;
};

#include "Mesh.ipp"

#endif


// dim      : 1D/2D/3D
// numNodes : number of nodes of the element
// beta     : cutoff ratio (element particle volume/element volume)
              //! which determines partially filled elements 
