/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         University Of Cambridge
File: Element.hpp
****************************************************************************/
#ifndef MPM_ELEMENT_H
#define MPM_ELEMENT_H

#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include <functional>

#include <eigen3/Eigen/Dense>

#include "Constants.hpp"
#include "Node.hpp"
//#include "element_matrix.hpp"

namespace mpm {
  class Element;
  template <class T>
  struct NewComparator {
    bool operator () (const T& arg1, const T& arg2) const {
      unsigned firstId = arg1 -> give_id();
      unsigned secondId = arg2 -> give_id();
      return firstId < secondId;
    }
  };
}

class mpm::Element {

protected:                    
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static const unsigned numMats = mpm::constants::NUMBODIES;
  static constexpr double beta = 0.7;
  static constexpr double alpha = 0.38;

public:
  // default constructor
  //! param[in] id element ID 
  //! param[in] mesh_space spacing of mesh taken equal to the element length 
  Element(const unsigned& id, const Eigen::Matrix<double,1,dim> &mesh_space);

  // initialise element
  void initialise_element();

  // insert the pointer to nodes of the element
  //! param[in] i index
  //! param[in] node_id id of the node
  //! param[in] node_ptr pointer to the node
  void set_element_nodes(const unsigned& i, const unsigned& node_id, mpm::Node* node_ptr) {
    elem_nodes_id_(i) = node_id;
    elem_nodes_ptrs_(i) = node_ptr;
  }

  // set the neighbouring elements
  //! param[in] elem_ptr pointer to the neighbour element
  void set_neighbour_elements(mpm::Element* elem_ptr) {
    neighbour_elements_.push_back(elem_ptr);
  }

  // compute centre coordinates
  void compute_centre_coordinates();

  // set the element status as free : NOT IN USE
  //! param[in] value value of the prescribed pressure at the free surface
  void set_element_nodes_pressure_constraints(const double& value) {
    for (unsigned i = 0; i < numNodes; i++)
      elem_nodes_ptrs_(i)->set_node_pressure_constraints(value);
    free_ = 1;
  }

  // set the number and volume of particles in the element
  //! param[in] p_volume volume of the particle
  void add_to_element_particle_density(const double& p_volume, const double& p_density) {
    elem_particle_density_ += 1;
    elem_material_volume_ += p_volume;
    elem_material_mass_ += (p_volume * p_density);
  }

  // find free surface nodes
  void find_free_surface_nodes(); 

  // give element id
  //! param[out] element id
  unsigned give_id() const {
    return elem_id_;
  }

  // give ids of element nodes
  //! param[out] vector of ids of element nodes
  Eigen::Matrix<unsigned,1,numNodes> give_element_node_ids() const {
    return elem_nodes_id_;
  }

  // give the vector of pointers to the nodes of the element
  //! param[out] vector of pointers to element nodes
  Eigen::Matrix<mpm::Node*, 1, numNodes> give_element_node_ptrs() const {
    return elem_nodes_ptrs_;
  }

  // give pointer to the node at given index
  //! param[out] pointer to a node at gien index
  mpm::Node* give_element_node_ptr_at(unsigned &index) const {
    return elem_nodes_ptrs_(index);
  }

  // give coordinates of the element centre
  //! param[out] element centre coordinates
  Eigen::Matrix<double,1,dim> give_element_centre_coord()const {
    return elem_centre_coord_;
  }

  // give element length
  //! param[out] element length
  Eigen::Matrix<double,1,dim> give_element_length() const {
    return elem_length_;
  }

  // give element volume
  //! param[out] element volume
  double give_element_volume() const {
    return elem_volume_;
  }

  // give pointers neighbour elements
  //! param[out] pointers to neighbour elements
  std::vector<mpm::Element*> give_neighbour_elements() {
    return neighbour_elements_;
  }

  // give particle density of the element
  //! param[out] number of particles in the element
  unsigned give_particle_density() const {
    return elem_particle_density_;
  }

  // give material volume in the element
  //! param[out] sum of volumes of particles in the element
  double give_material_volume_in_element() const {
    return elem_material_volume_;
  }

  // give material volume / element volume ratio
  //! param[out] ratio volume of particles in element / element volume
  double give_ratio_particle_volume_to_element_volume() const {
    double ratio = elem_material_volume_ / elem_volume_;
    return ratio;
  }

  double give_ratio_particle_mass_to_element_mass() const {
    double ratio = elem_material_mass_ / (elem_volume_ * 2600);
    return ratio;
  }

  // give the status of partially/fully filled with particles
  //! param[out] status true/false 
  bool is_element_fully_filled() {
    if ((elem_material_volume_ / elem_volume_) < beta)
      return 0;
    else
      return 1;
  }

  // give status of free element
  //! param[out] status of free
  bool give_status_of_free_element() const {
    return free_;
  }

  // give pointers to nodes located at the free surface
  void give_element_free_surface_nodes(std::vector<mpm::Node*> &free_nodes) {
    if (free_) {
      for (const auto &n_ptr : free_surface_nodes_)
	free_nodes.push_back(n_ptr);
    }
  }

  // build element mass matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value_ss entry for solid at (row,column) 
  //! param[in] value_ff entry for water at (row,column) 
  void build_M_matrix_mp(const unsigned &row, const unsigned& column, const double& value_ss, const double& value_ff);

  // build element K bar matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_bar_matrix_mp(const unsigned &row, const unsigned& column, const double& value);

  // build element Laplace matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_Laplace_matrix_mp(const unsigned &row, const unsigned& column, const double& value);

  // build element K_L2 matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_L2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value);

  // build element K_L3 matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_L3_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value);

  // build element K_L4 matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_L4_matrix_mp(const unsigned& row, const unsigned& column, const double &value);

  // build element K_S2 matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_S2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value);

  // build element K_W2 matrix evaluated at material points
  //! param[in] row
  //! param[in] column
  //! param[in] value entry at (row,column) 
  void build_K_W2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value);

  // compute element matrix evaluated at gauss points
  void compute_element_matrix_gp(const double& permeability);


  // give solid mass matrix
  //! param[out] M_SS element solid mass matrix
  Eigen::Matrix<double,numNodes,numNodes> give_M_SS_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta) 
    return M_SS_mp_; 
    // else
    // return M_SS_gp_;
  }

  // give water mass matrix
  //! param[out] M_FF element water mass matrix
  Eigen::Matrix<double,numNodes,numNodes> give_M_FF_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta) 
    return M_FF_mp_; 
    // else
    // return M_FF_gp_;
  }

  // give K bar matrix
  //! param[out] K_BAR element K bar matrix
  Eigen::Matrix<double,numNodes,numNodes> give_K_BAR_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
    return K_BAR_mp_;
    // else
    // return K_BAR_gp_;
  }

  // give Laplace matrix
  //! param[out] Laplace_ element laplace matrix
  Eigen::Matrix<double,numNodes,numNodes> give_L_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
      return L_mp_;
      // else
      // return L_gp_;
  }

  // give K_L2 matrix
  //! param[out] K_L2 element K_L2 matrix
  Eigen::Matrix<double,numNodes,dim*numNodes> give_K_L2_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
    return K_L2_mp_;
    // else
    //   return K_L2_gp_;
  }

  // give K_L3 matrix
  //! param[out] K_L3 element K_L2 matrix
  Eigen::Matrix<double,numNodes,dim*numNodes> give_K_L3_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
      return K_L3_mp_;
      // else
      // return K_L3_gp_;
  }

  // give K_L4 matrix
  //! param[out] K_L4 element K_L4 matrix
  Eigen::Matrix<double,numNodes,numNodes> give_K_L4_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
      return K_L4_mp_;
      // else
      // return K_L4_gp_;
  }

  // give K_S2 matrix
  //! param[out] K_S2_
  Eigen::Matrix<double,dim*numNodes,numNodes> give_K_S2_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
      return K_S2_mp_;
      // else
      // return K_S2_gp_;
  }

  // give K_W2 matrix
  //! param[out] K_W2_
  Eigen::Matrix<double,dim*numNodes,numNodes> give_K_W2_matrix() const {
    //if (std::fabs(elem_material_volume_ / elem_volume_) < beta)
      return K_W2_mp_;
      // else
      // return K_W2_gp_;
  }

  void element_restrained() {
    restrained_ = true;
  }

private:
    // compare the nodes of two elements
    // and insert the common nodes to vector of free_nodes
    //! param[in] elem_ptr pointer to the element which is compared
    void insert_to_free_surface_nodes(const mpm::Element* elem_ptr);

protected:
  // element ID
  unsigned elem_id_;

  // IDs of the nodes of the element
  Eigen::Matrix<unsigned, 1, numNodes> elem_nodes_id_;
  // pointers to the nodes of the element
  Eigen::Matrix<mpm::Node*, 1, numNodes> elem_nodes_ptrs_;
  // pointers to the neighbour elements
  std::vector<mpm::Element*> neighbour_elements_;

  // Element length
  Eigen::Matrix<double, 1, dim> elem_length_;
  // Element centre coordinates
  Eigen::Matrix<double, 1, dim> elem_centre_coord_;

  // Element Volume 
  double elem_volume_;
  // Particle Density
  unsigned elem_particle_density_;
  // Material volume : sum of volumes of particles in the element
  double elem_material_volume_;
  double elem_material_mass_;
  // status of free element or not 
  bool free_;
  bool restrained_;
  // pointers to the free surface nodes of the element
  std::set<mpm::Node*, mpm::NewComparator<mpm::Node*> > free_surface_nodes_;

  // shared pointer to the element matrix
  //std::shared_ptr<mpm::ElementMatrix> element_matrix_ptr_;

  // element mass matrix
  Eigen::Matrix<double, numNodes, numNodes> M_SS_mp_;
  Eigen::Matrix<double, numNodes, numNodes> M_FF_mp_;
  Eigen::Matrix<double, numNodes, numNodes> M_SS_gp_;
  Eigen::Matrix<double, numNodes, numNodes> M_FF_gp_;
  // element K_bar matrix
  Eigen::Matrix<double, numNodes, numNodes> K_BAR_mp_;
  Eigen::Matrix<double, numNodes, numNodes> K_BAR_gp_;
  std::array<Eigen::Matrix<double, numNodes, numNodes>, dim> gauss_D_;

  Eigen::Matrix<double, numNodes, numNodes> L_gp_;
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L2_gp_;
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L3_gp_;
  // ******************************************************
  // element Laplace matrix
  Eigen::Matrix<double, numNodes, numNodes> L_mp_;
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L2_mp_;
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L3_mp_;
  Eigen::Matrix<double, dim*numNodes, numNodes> K_S2_mp_;
  // ******************************************************
  
  Eigen::Matrix<double, numNodes, numNodes> K_L4_mp_;
  Eigen::Matrix<double, numNodes, numNodes> K_L4_gp_;

  Eigen::Matrix<double, dim*numNodes, numNodes> K_S2_gp_;
  Eigen::Matrix<double, dim*numNodes, numNodes> K_W2_mp_;
  Eigen::Matrix<double, dim*numNodes, numNodes> K_W2_gp_;
};

#include "Element.ipp"

#endif

// dim      : 1D/2D/3D
// numNodes : number of nodes of the element
// beta     : cutoff ratio (element particle volume/element volume)
              //! which determines the element partially/fully filled 
// alpha    : cutoff ratio (element particle volume/element volume)
              //! which determines free surface nodes
