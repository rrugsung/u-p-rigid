 static const unsigned dim = mpm::constants::DIM;
 static const unsigned numNodes = mpm::constants::NUMNODES;
mpm::ElementMatrix::ElementMatrix(const Eigen::Matrix<double,1,dim> element_length) {
    length_ = element_length;
    const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
    Eigen::Vector2d qpoint_one(-one_by_sqrt3, -one_by_sqrt3);
    Eigen::Vector2d qpoint_two(one_by_sqrt3, -one_by_sqrt3);
    Eigen::Vector2d qpoint_three(one_by_sqrt3, one_by_sqrt3);
    Eigen::Vector2d qpoint_four(-one_by_sqrt3, one_by_sqrt3);

    sfun_.at(0) = shape_fun(qpoint_one);
    sfun_.at(1) = shape_fun(qpoint_two);
    sfun_.at(2) = shape_fun(qpoint_three);
    sfun_.at(3) = shape_fun(qpoint_four);

    grad_sfun_.at(0) = grad_shape_fun(qpoint_one);
    grad_sfun_.at(1) = grad_shape_fun(qpoint_two);
    grad_sfun_.at(2) = grad_shape_fun(qpoint_three);
    grad_sfun_.at(3) = grad_shape_fun(qpoint_four);
}

Eigen::Matrix<double,1,numNodes> mpm::ElementMatrix::shape_fun(const Eigen::Vector2d &xi) {
    Eigen::Matrix<double,1,numNodes> s_fun;
    s_fun(0) = 0.25 * std::fabs((1 - xi(0)) * (1 - xi(1)));
    s_fun(1) = 0.25 * std::fabs((1 + xi(0)) * (1 - xi(1)));
    s_fun(2) = 0.25 * std::fabs((1 + xi(0)) * (1 + xi(1)));
    s_fun(3) = 0.25 * std::fabs((1 - xi(0)) * (1 + xi(1)));
    return s_fun;
}


Eigen::Matrix<double,dim,numNodes> mpm::ElementMatrix::grad_shape_fun(const Eigen::Vector2d &xi) {
    Eigen::Matrix<double, dim, numNodes> grad_s_fun;
    grad_s_fun(0, 0) = -0.5 * (1 - xi(1)) / length_(0);
    grad_s_fun(0, 1) =  0.5 * (1 - xi(1)) / length_(0);
    grad_s_fun(0, 2) =  0.5 * (1 + xi(1)) / length_(0);
    grad_s_fun(0, 3) = -0.5 * (1 + xi(1)) / length_(0);

    grad_s_fun(1, 0) = -0.5 * (1 - xi(0)) / length_(1);
    grad_s_fun(1, 1) = -0.5 * (1 + xi(0)) / length_(1);
    grad_s_fun(1, 2) =  0.5 * (1 + xi(0)) / length_(1);
    grad_s_fun(1, 3) =  0.5 * (1 - xi(0)) / length_(1);
    return grad_s_fun;
}

void mpm::ElementMatrix::assign_porosity_to_gauss_points(const Eigen::Matrix<double,1,numNodes>& nodal_solid_porosity, const Eigen::Matrix<double,1,numNodes>& nodal_water_porosity) {
  solid_porosity_gp = Eigen::Matrix<double,1,4>::Zero();
  water_porosity_gp = Eigen::Matrix<double,1,4>::Zero();
  for(unsigned i = 0; i < 4; i++) {
    Eigen::Matrix<double, 1, numNodes> sfun_gp = sfun_.at(i);
    for (unsigned j = 0; j < numNodes; j++) {
      solid_porosity_gp(i) += sfun_gp(j) * nodal_solid_porosity(j);
      water_porosity_gp(i) += sfun_gp(j) * nodal_water_porosity(j);
    }
  }
}

Eigen::Matrix<double, numNodes, numNodes> mpm::ElementMatrix::M_ss_matrix_gp() {
  M_ss_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double, 1, numNodes> qp_sfun = sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        M_ss_(i,j) += (1.0 * solid_porosity_gp(q) * solidgraindensity * qp_sfun(i) * qp_sfun(j));
      }
    }
  }
  return M_ss_;
}

Eigen::Matrix<double, numNodes, numNodes> mpm::ElementMatrix::M_ff_matrix_gp() {
  M_ff_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double, 1, numNodes> qp_sfun = sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        M_ff_(i,j) += (1.0 * water_porosity_gp(q) * watergraindensity * qp_sfun(i) * qp_sfun(j));
      }
    }
  }
  return M_ff_;
}


Eigen::Matrix<double, numNodes, numNodes> mpm::ElementMatrix::K_bar_matrix_gp(const double permeability) {
  Kbar_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        Kbar_(i,j) += 1.0 * ((water_porosity_gp(q) * water_porosity_gp(q)) * watergraindensity * 9.81 / permeability) * qp_sfun(i) * qp_sfun(j);
        // Kbar_(i,j) += 1.0 * ((water_porosity_gp(q) * water_porosity_gp(q)) * watergraindensity * 9.81) * qp_sfun(i) * qp_sfun(j);
      }
    }
  }
  // std::cout << Kbar_ << "\n";
  return Kbar_;
}

Eigen::Matrix<double, numNodes, numNodes> mpm::ElementMatrix::Laplace_matrix_gp() {
  Laplace_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double, dim, numNodes> qp_grad_sfun = grad_sfun_.at(q);
    double laplace_constant = ((solid_porosity_gp(q) / solidgraindensity) + (water_porosity_gp(q) / watergraindensity));
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        Laplace_(i,j) += 1.0 * laplace_constant * ((qp_grad_sfun(0,i) * qp_grad_sfun(0,j)) + (qp_grad_sfun(1,i) * qp_grad_sfun(1,j)));
      }
    }
  }
  return Laplace_;
}

Eigen::Matrix<double, numNodes, dim*numNodes> mpm::ElementMatrix::K_L2_matrix_gp() {
  K_L2_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    Eigen::Matrix<double, dim, numNodes> qp_grad_sfun = grad_sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        K_L2_(i,j) += 1.0 * water_porosity_gp(q) * qp_grad_sfun(0,i) * qp_sfun(j);
        K_L2_(i, (j + numNodes)) += 1.0 * water_porosity_gp(q) * qp_grad_sfun(1,i) * qp_sfun(j);
      }
    }
  }
  return K_L2_;
}

Eigen::Matrix<double, numNodes, dim*numNodes> mpm::ElementMatrix::K_L3_matrix_gp() {
  K_L3_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    Eigen::Matrix<double, dim, numNodes> qp_grad_sfun = grad_sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        K_L3_(i,j) += -1.0 * qp_sfun(i) * qp_grad_sfun(0,j) ;
        K_L3_(i, (j + numNodes)) += -1.0 * qp_sfun(i) * qp_grad_sfun(1,j) ;
      }
    }
  }
  return K_L3_;
}

Eigen::Matrix<double, numNodes, numNodes> mpm::ElementMatrix::K_L4_matrix_gp() {
  K_L4_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        K_L4_(i,j) += qp_sfun(i) * qp_sfun(j);
      }
    }
  }
  return K_L4_;
}

Eigen::Matrix<double, dim*numNodes, numNodes> mpm::ElementMatrix::K_S2_matrix_gp() {
  K_S2_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    Eigen::Matrix<double, dim, numNodes> qp_grad_sfun = grad_sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        K_S2_(i,j) += -1.0 * solid_porosity_gp(q) * qp_sfun(i) * qp_grad_sfun(0,j) ;
        K_S2_((i + numNodes), j) += -1.0 * solid_porosity_gp(q) * qp_sfun(i) * qp_grad_sfun(1,j) ;
        // K_S2_(i,j) += -1.0 * qp_sfun(i) * qp_grad_sfun(0,j) ;
        // K_S2_((i + numNodes), j) += -1.0 * qp_sfun(i) * qp_grad_sfun(1,j) ;
      }
    }
  }
  return K_S2_;
}

Eigen::Matrix<double, dim*numNodes, numNodes> mpm::ElementMatrix::K_W2_matrix_gp() {
  K_W2_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  for (unsigned q = 0; q < 4; q++) {
    Eigen::Matrix<double,1,numNodes> qp_sfun = sfun_.at(q);
    Eigen::Matrix<double, dim, numNodes> qp_grad_sfun = grad_sfun_.at(q);
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        K_W2_(i,j) += -1.0 * water_porosity_gp(q) * qp_sfun(i) * qp_grad_sfun(0,j) ;
        K_W2_((i + numNodes), j) += -1.0 * water_porosity_gp(q) * qp_sfun(i) * qp_grad_sfun(1,j);
      }
    }
  }
  return K_W2_;
}


void mpm::ElementMatrix::initialise_element_matrix() {
  M_ss_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_ff_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  Kbar_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  Laplace_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_L2_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L3_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_S2_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_W2_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
}



