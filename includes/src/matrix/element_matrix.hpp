/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         University Of Cambridge
File: ElementMatrix.hpp
****************************************************************************/
#ifndef MPM_ELEMENTMATRIX_H
#define MPM_ELEMENTMATRIX_H
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "Constants.hpp"

namespace mpm {
    class ElementMatrix;
}

// Quadrilateral Element Matrix Class
//! stores element matrices
//! computed at Gaussian quadrature points
class mpm::ElementMatrix {
protected:
  static const unsigned dim = mpm::constants::DIM;
  static const unsigned numNodes = mpm::constants::NUMNODES;
  static constexpr double solidgraindensity = 2600;
  static constexpr double watergraindensity = 1000;
  //static constexpr double permeability = 0.00000001;

public:

    // constructor
    //! initialise the shape functions and gradient shape functions
    //! param[in] element_length 
    ElementMatrix(const Eigen::Matrix<double,1,dim> element_length);

    // initialise matrices
    void initialise_element_matrix();

    // shape functions
    //! compute and return shape functions
    //! param[in] xi local coordinates of the quadrature point
    //! param[out] shape functions at given quadrature point
    Eigen::Matrix<double,1,numNodes> shape_fun(const Eigen::Vector2d &xi);

    // global gradient of shape functions
    //! compute and return global gradient of shape functions
    //! param[in] xi local coordinates of the quadrature point
    //! param[out] gradient of shape functions at given quadrature point
    Eigen::Matrix<double,dim,numNodes> grad_shape_fun(const Eigen::Vector2d &xi);

  // assign solid and water porosity values to gauss points
  //! param[in] nodal solid poroisty
  //! param[in] nodal water porosity
  void assign_porosity_to_gauss_points(const Eigen::Matrix<double,1,numNodes>& nodal_solid_porosity, const Eigen::Matrix<double,1,numNodes>& nodal_water_porosity);

  // evaluate the mass matrix
  //! param[out] M_
  Eigen::Matrix<double, numNodes, numNodes> M_ss_matrix_gp();

  // evaluate laplace matrix
  //! param[out] L_
  Eigen::Matrix<double, numNodes, numNodes> M_ff_matrix_gp();

  // evaluate divergence/gradient matrix
  //! param[out] D_
  Eigen::Matrix<double, numNodes, numNodes> K_bar_matrix_gp(const double permeability);

  // evaluate laplace matrix
  //! param[out] Laplace_
  Eigen::Matrix<double, numNodes, numNodes> Laplace_matrix_gp();

  // evaluate K_L2 matrix
  //! param[out] K_L2_
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L2_matrix_gp();

  // evaluate K_L3 matrix
  //! param[out] K_L3_
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L3_matrix_gp();

  // evaluate K_L4 matrix
  //! param[out] K_L3_
  Eigen::Matrix<double, numNodes, numNodes> K_L4_matrix_gp();

  // evaluate K_S2 matrix
  //! param[out] K_S2_
  Eigen::Matrix<double, dim*numNodes, numNodes> K_S2_matrix_gp();

  // evaluate K_W2 matrix
  //! param[out] K_W2_
  Eigen::Matrix<double, dim*numNodes, numNodes> K_W2_matrix_gp();



protected:
  // element length
  Eigen::Matrix<double,1,dim> length_;
  // shape functions at quadrature points
  // shape functions at quadrature points
  std::array<Eigen::Matrix<double, 1, numNodes>,4> sfun_;
  // global gradient of shape functions
  std::array<Eigen::Matrix<double, dim, numNodes>,4> grad_sfun_;

  // solid porosity
  Eigen::Matrix<double,1,4> solid_porosity_gp;
  Eigen::Matrix<double,1,4> water_porosity_gp;

  // element mass matrix computed at gauss points
  Eigen::Matrix<double, numNodes, numNodes> M_ss_;
  Eigen::Matrix<double, numNodes, numNodes> M_ff_;
  // element K bar matrix computed at gauss points
  Eigen::Matrix<double, numNodes, numNodes> Kbar_;

  // element laplace matrix computed at gauss points
  Eigen::Matrix<double, numNodes, numNodes> Laplace_;

  Eigen::Matrix<double, numNodes, dim*numNodes> K_L2_;
  Eigen::Matrix<double, numNodes, dim*numNodes> K_L3_;
  Eigen::Matrix<double, numNodes, numNodes> K_L4_;

  Eigen::Matrix<double, dim*numNodes, numNodes> K_S2_;
  Eigen::Matrix<double, dim*numNodes, numNodes> K_W2_;
};

#include "element_matrix.ipp"

#endif
