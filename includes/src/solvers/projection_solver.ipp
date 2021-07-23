
mpm::ProjectionSolver::ProjectionSolver(const double &permeability, const double &compressibility) {
    element_ptrs_.clear();
    node_ptrs_.clear();
    velocity_dof_ = 0;
    pressure_dof_ = 0;
    const_permeability_ = permeability;
    fluid_compressibility_ = compressibility;
    std::cout << const_permeability_ << "\n";
    std::cout << fluid_compressibility_ << "\n";
}

void mpm::ProjectionSolver::assemble_solver(mpm::Mesh* mesh_ptr) {
    element_ptrs_.clear();
    node_ptrs_.clear();
    std::set<mpm::Node*, mpm::setcomparator<mpm::Node*> > node_ptrs_set;
    node_ptrs_set.clear();
    double cut_off = theta;
    mesh_ptr->give_elements_with_particles_above_cutoff(element_ptrs_, cut_off);
    for (const auto &elem_ptr : element_ptrs_) {
        Eigen::Matrix<mpm::Node*, 1, numNodes> elem_nodes = elem_ptr->give_element_node_ptrs();
        for (unsigned i = 0; i < numNodes; i++)
          node_ptrs_set.insert(elem_nodes(i));
    }
    for (const auto &node_ptr : node_ptrs_set){
        node_ptrs_.push_back(node_ptr);
    }

    num_nodes_ = node_ptrs_.size();
    velocity_dof_ = num_nodes_ * dim;
    pressure_dof_ = num_nodes_;
    this->element_mapping_matrix();
    this->store_pressure_constraints();

    L_.resize(pressure_dof_,pressure_dof_);
    L_.setZero();
    L_.reserve(Eigen::VectorXd::Constant(num_nodes_, 8));
    this->build_L_matrix();

    K_L2_.resize(num_nodes_, velocity_dof_);
    K_L2_.setZero();
    K_L2_.reserve(Eigen::VectorXd::Constant(velocity_dof_, 16));
    K_L3_.resize(num_nodes_, velocity_dof_);
    K_L3_.setZero();
    K_L3_.reserve(Eigen::VectorXd::Constant(velocity_dof_, 16));
    this->build_K_L_matrix();

    K_S2_.resize(velocity_dof_, pressure_dof_);
    K_S2_.setZero();
    K_S2_.reserve(Eigen::VectorXd::Constant(pressure_dof_, 8));
    this->build_K_S2_matrix();
}

bool mpm::ProjectionSolver::solve_pressure_poisson_equation(const double& dt) {
  stiffness_matrix_ =  L_;

  intermediate_solid_velocity_.resize(velocity_dof_);
  intermediate_water_velocity_.resize(velocity_dof_);
  unsigned node = 0;
  
  for(const auto nptr : node_ptrs_) {
    Eigen::Matrix<double,1,dim> nsolid_int_velocity = nptr->give_node_solid_intermediate_velocity();
    intermediate_solid_velocity_(node) = nsolid_int_velocity(0);
    intermediate_solid_velocity_(node+num_nodes_) = nsolid_int_velocity(1);
    Eigen::Matrix<double,1,dim> nwater_int_velocity = nptr->give_node_multimaterial_water_intermediate_velocities(0);
    intermediate_water_velocity_(node) = nwater_int_velocity(0);
    intermediate_water_velocity_(node+num_nodes_) = nwater_int_velocity(1);
    node++;
  }

  force_vector_ = (K_L2_ * intermediate_solid_velocity_) - (K_L3_ * intermediate_water_velocity_); 
  this->apply_pressure_boundary_conditions_to_system();
  if (!Conjugate_Gradient(pressure_))
    std::cerr << "Failed solving Poisson equation using CG" << "\n";
  std::cout << "force_vector (399): \n" << force_vector_(399) << "\n";
  std::cout << "force_vector (400): \n" << force_vector_(400) << "\n";
  std::cout << "force_vector (420): \n" << force_vector_(420) << "\n";
  std::cout << "force_vector (421): \n" << force_vector_(421) << "\n";
  //std::cout << "L_: \n" << stiffness_matrix_ << "\n";
  std::cout << "pressure_ (399): \n" << pressure_(399) <<  "\n";
  std::cout << "pressure_ (400): \n" << pressure_(400) <<  "\n";
  std::cout << "pressure_ (420): \n" << pressure_(420) <<  "\n";
  std::cout << "pressure_ (421): \n" << pressure_(421) <<  "\n";
  //for (unsigned i = 0; i < pressure_.size(); i++) {
  //   if (std::fabs(pressure_(i)) < 1.0E2)
  //     pressure_(i) = 0.;
  // }
}


bool mpm::ProjectionSolver::solve_final_acceleration(const double &dt, const double &gamma) {

}

void mpm::ProjectionSolver::element_mapping_matrix() {
  unsigned num_elem = element_ptrs_.size();
  element_map_.resize(num_elem, numNodes);

  std::vector<mpm::Node*>::iterator it;
  unsigned element = 0;
  for (const auto &elem_ptr : element_ptrs_) {
    Eigen::Matrix<mpm::Node*,1,numNodes> elem_nodes = elem_ptr->give_element_node_ptrs();
    for (unsigned j = 0; j < numNodes; j++) {
      it = std::find(node_ptrs_.begin(), node_ptrs_.end(), elem_nodes(j));
      if(it != node_ptrs_.end())
	element_map_(element,j) = it - node_ptrs_.begin();
    }
    element++;
  }
}


void mpm::ProjectionSolver::store_pressure_constraints() {
  pressure_constraints_.clear();
  for (const auto &node_ptr : node_ptrs_) {
    auto constraint = node_ptr->give_node_pressure_constraint();
    if (mpm::misc::scalar_beta_)
      std::get<1>(constraint) = 0.;
    pressure_constraints_.push_back(constraint);
  }
  if (pressure_constraints_.size() != pressure_dof_) {
    std::cerr << "Error file: projection_solver.ipp Line: 60" << "\n";
    std::cerr << "\t sizes of pressure constraints and pressure degrees of freedom are mismatched" << "\n";
    abort();
  }
}

void mpm::ProjectionSolver::build_M_matrix() {
    unsigned elem = 0;
    Eigen::Matrix<double,numNodes,numNodes> elem_M_11;
    Eigen::Matrix<double,numNodes,numNodes> elem_M_22;
    for (const auto &elem_ptr : element_ptrs_) {
        elem_M_11 = elem_ptr->give_M_SS_matrix();
        elem_M_22 = elem_ptr->give_M_FF_matrix();
        for (unsigned i = 0; i < numNodes; i++) {
            for (unsigned j = 0; j < numNodes; j++) {
                for (unsigned k = 0; k < dim; k++) {
                    unsigned row = element_map_(elem, i) + (num_nodes_ * k);
                    unsigned col = element_map_(elem, j) + (num_nodes_ * k);
                    M_11_.coeffRef(row,col) += elem_M_11(i,j);
                    M_22_.coeffRef(row,col) += elem_M_22(i,j);
                    M_.coeffRef(row,col) += elem_M_11(i,j);
                    M_.coeffRef(row,(col + velocity_dof_)) += elem_M_22(i,j);
                    M_.coeffRef((row + velocity_dof_), (col + velocity_dof_)) += elem_M_22(i,j);
                }
            }
        }
        elem++;
    }
}


void mpm::ProjectionSolver::build_Kbar_matrix() {
    unsigned elem = 0;
    Eigen::Matrix<double,numNodes,numNodes> elem_K_bar;
    for (const auto &elem_ptr : element_ptrs_) {
        elem_K_bar = elem_ptr->give_K_BAR_matrix();
        for (unsigned i = 0; i < numNodes; i++) {
            for (unsigned j = 0; j < numNodes; j++) {
                for (unsigned k = 0; k < dim; k++) {
                    unsigned row = element_map_(elem, i) + (num_nodes_ * k);
                    unsigned col = element_map_(elem, j) + (num_nodes_ * k);
                    // K_bar_.coeffRef(row,col) += elem_K_bar(i,j);
                    // K_bar_.coeffRef(row, (col + velocity_dof_)) += -elem_K_bar(i,j);
                    K_bar_.coeffRef((row + velocity_dof_), col) += -elem_K_bar(i,j);
                    K_bar_.coeffRef((row + velocity_dof_), (col + velocity_dof_)) += elem_K_bar(i,j);
                    K_single_.coeffRef(row,col) += elem_K_bar(i,j);
                }
            }
        }
        elem++;
    }
}


void mpm::ProjectionSolver::build_L_matrix() {
  unsigned elem = 0;
  Eigen::Matrix<double,numNodes,numNodes> elem_L;
  for (const auto &elem_ptr : element_ptrs_) {
    elem_L = elem_ptr->give_L_matrix();
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        unsigned row = element_map_(elem, i);
        unsigned col = element_map_(elem, j);
        L_.coeffRef(row,col) += elem_L(i,j);
      }
    }
    elem++;
  }
}

void mpm::ProjectionSolver::build_K_L_matrix() {
  unsigned elem = 0;
  Eigen::Matrix<double,numNodes,dim*numNodes> elem_K_L2, elem_K_L3;
  for (const auto &elem_ptr : element_ptrs_) {
    elem_K_L2 = elem_ptr->give_K_L2_matrix();
    elem_K_L3 = elem_ptr->give_K_L3_matrix();
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        for (unsigned k = 0; k < dim; k++) {
          unsigned row = element_map_(elem, i);
          unsigned col = element_map_(elem, j) + (num_nodes_ * k);
          K_L2_.coeffRef(row,col) += elem_K_L2(i,(j + k*numNodes));
          K_L3_.coeffRef(row,col) += elem_K_L3(i,(j + k*numNodes));
        }
      }
    }
    elem++;
  }
}

void mpm::ProjectionSolver::build_K_L4_matrix() {
  unsigned elem = 0;
  Eigen::Matrix<double,numNodes,numNodes> elem_K_L4;
  for (const auto &elem_ptr : element_ptrs_) {
    elem_K_L4 = elem_ptr->give_K_L4_matrix();
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        unsigned row = element_map_(elem, i);
        unsigned col = element_map_(elem, j);
        K_L4_.coeffRef(row,col) += elem_K_L4(i,j);
      }
    }
    elem++;
  }
}

void mpm::ProjectionSolver::build_K_S2_matrix() {
  unsigned elem = 0;
  Eigen::Matrix<double,dim*numNodes,numNodes> elem_K_S2;
  for (const auto &elem_ptr : element_ptrs_) {
    elem_K_S2 = elem_ptr->give_K_S2_matrix();
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        for (unsigned k = 0; k < dim; k++) {
          unsigned row = element_map_(elem, i) + (num_nodes_ * k);
          unsigned col = element_map_(elem, j);
          K_S2_.coeffRef(row,col) += elem_K_S2((i + k*numNodes),j);
        }
      }
    }
    elem++;
  }
}

void mpm::ProjectionSolver::build_K_W2_matrix() {
  unsigned elem = 0;
  Eigen::Matrix<double,dim*numNodes,numNodes> elem_K_W2;
  for (const auto &elem_ptr : element_ptrs_) {
    elem_K_W2 = elem_ptr->give_K_W2_matrix();
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j = 0; j < numNodes; j++) {
        for (unsigned k = 0; k < dim; k++) {
          unsigned row = element_map_(elem, i) + (num_nodes_ * k);
          unsigned col = element_map_(elem, j);
          K_W2_.coeffRef(row,col) += elem_K_W2((i + k*numNodes),j);
        }
      }
    }
    elem++;
  }
}


bool mpm::ProjectionSolver::Conjugate_Gradient(Eigen::VectorXd &solution) {
  this->check_for_zeroes_in_linear_system();
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> CG_solver;
  CG_solver.setMaxIterations(100*stiffness_matrix_.cols());
  CG_solver.setTolerance(0.000001);
  solution = CG_solver.compute(stiffness_matrix_).solve(force_vector_);
  if (CG_solver.info() != Eigen::Success)
    return 0;
  else
    return 1;
}

bool mpm::ProjectionSolver::Least_Squares_Conjugate_Gradient(Eigen::VectorXd &solution) {
    this->check_for_zeroes_in_linear_system();
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> LSCG_solver;
    LSCG_solver.setMaxIterations(100*stiffness_matrix_.cols());
    LSCG_solver.setTolerance(0.000001);
    solution = LSCG_solver.compute(stiffness_matrix_).solve(force_vector_);
    if (LSCG_solver.info() != Eigen::Success)
        return 0;
    else
        return 1;
}

void mpm::ProjectionSolver::apply_pressure_boundary_conditions_to_system() {
    if((pressure_constraints_.size() != stiffness_matrix_.rows()) || (pressure_constraints_.size() != stiffness_matrix_.cols())) {
        std::cerr << "Error file: projection_solver.ipp Line: 163" << "\n";
        std::cerr << "\t sizes of stifness matrix and pressure constraints are mismatched" << "\n";
        abort();
    }
    if (pressure_constraints_.size() != force_vector_.size()) {
        std::cerr << "Error file: projection_solver.ipp Line: 168" << "\n";
        std::cerr << "\t sizes of force vector and pressure constraints are mismatched" << "\n";
        abort();
    }
    unsigned n = 0;
    for (const auto &npressure_bc : pressure_constraints_) {
        if(std::get<0>(npressure_bc)) {
            for (unsigned p = 0; p < force_vector_.size(); p++) 
              force_vector_(p) -= (std::get<1>(npressure_bc)) * stiffness_matrix_.coeffRef(p,n);
            // force_vector_(p) -= 0.0 * stiffness_matrix_.coeffRef(p,n);
            for (unsigned i = 0; i < pressure_dof_; i++) {
                stiffness_matrix_.coeffRef(n,i) = 0;
                stiffness_matrix_.coeffRef(i,n) = 0;
            }
            force_vector_(n) = std::get<1>(npressure_bc);
            if (mpm::misc::scalar_beta_)
              force_vector_(n) = 0.;
            // force_vector_(n) = 0.0;
            stiffness_matrix_.coeffRef(n,n)= 1;
        }   
        n++;
    }
}


void mpm::ProjectionSolver::apply_velocity_boundary_conditions_to_system() {
    if((velocity_constraints_.size() != stiffness_matrix_.rows()) || (velocity_constraints_.size() != stiffness_matrix_.cols())) {
        std::cerr << "Error file: projection_solver.ipp Line: 348" << "\n";
        std::cerr << "\t sizes of stifness matrix and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }
    if (velocity_constraints_.size() != force_vector_.size()) {
        std::cerr << "Error file: projection_solver.ipp Line: 353" << "\n";
        std::cerr << "\t sizes of force vector and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }

 typedef Eigen::Triplet<double> T;
  std::vector<T> update;
  update.reserve(stiffness_matrix_.nonZeros() +
                 velocity_constraints_.size() * (2 * velocity_dof_ + 1));

  std::map<std::pair<unsigned, unsigned>, bool> update_vector;

  unsigned n = 0;
  for (const auto &nvelocity_bc : velocity_constraints_) {
    if (std::get<0>(nvelocity_bc)) {
      const auto velocity_bc = std::get<1>(nvelocity_bc);
      force_vector_ -=  velocity_bc * stiffness_matrix_.col(n);
      force_vector_[n] = velocity_bc;

      for (unsigned i = 0; i < 2 * velocity_dof_; ++i) {
        update.push_back(T(n, i, 0));
        update.push_back(T(i, n, 0));
        update_vector.emplace(std::make_pair(std::make_pair(n, i), false));
        update_vector.emplace(std::make_pair(std::make_pair(i, n), false));
        // stiffness_matrix_.coeffRef(n, i) = 0;
        // stiffness_matrix_.coeffRef(i, n) = 0;
      }
      // stiffness_matrix_.coeffRef(n, n) = 1;
      update.push_back(T(n, n, 1));
      update_vector.emplace(std::make_pair(std::make_pair(n, n), false));
    }
    ++n;
  }

  // std::sort(update_vector.begin(), update_vector.end());
  for (int k = 0; k < stiffness_matrix_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness_matrix_, k);
         it; ++it) {
      if (update_vector.find(std::make_pair(it.row(), it.col())) ==
          update_vector.end())
        update.push_back(T(it.row(), it.col(), it.value()));
    }
  }
  stiffness_matrix_.setFromTriplets(update.begin(), update.end());
  stiffness_matrix_.makeCompressed();
}

void mpm::ProjectionSolver::apply_solid_velocity_boundary_conditions_to_system() {
    if((solid_velocity_constraints_.size() != stiffness_matrix_.rows()) || (solid_velocity_constraints_.size() != stiffness_matrix_.cols())) {
        std::cerr << "Error file: projection_solver.ipp Line: 376" << "\n";
        std::cerr << "\t sizes of stifness matrix and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }
    if (solid_velocity_constraints_.size() != force_vector_.size()) {
        std::cerr << "Error file: projection_solver.ipp Line: 381" << "\n";
        std::cerr << "\t sizes of force vector and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }

  typedef Eigen::Triplet<double> T;
  std::vector<T> update;
  update.reserve(stiffness_matrix_.nonZeros() +
                 velocity_constraints_.size() * (2 * velocity_dof_ + 1));

  std::map<std::pair<unsigned, unsigned>, bool> update_vector;

  unsigned n = 0;
  for (const auto &nvelocity_bc : solid_velocity_constraints_) {
    if (std::get<0>(nvelocity_bc)) {
      const auto velocity_bc = std::get<1>(nvelocity_bc);
      force_vector_ -=  velocity_bc * stiffness_matrix_.col(n);
      force_vector_[n] = velocity_bc;
      for (unsigned i = 0; i < velocity_dof_; i++) {
        // stiffness_matrix_.coeffRef(n, i) = 0;
        // stiffness_matrix_.coeffRef(i, n) = 0;
        update.push_back(T(n, i, 0));
        update.push_back(T(i, n, 0));
        update_vector.emplace(std::make_pair(std::make_pair(n, i), false));
        update_vector.emplace(std::make_pair(std::make_pair(i, n), false));
      }
      update.push_back(T(n, n, 1));
      update_vector.emplace(std::make_pair(std::make_pair(n, n), false));
    }
    n++;
  }
  // std::sort(update_vector.begin(), update_vector.end());
  for (int k = 0; k < stiffness_matrix_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness_matrix_, k);
         it; ++it) {
      if (update_vector.find(std::make_pair(it.row(), it.col())) ==
          update_vector.end())
        update.push_back(T(it.row(), it.col(), it.value()));
    }
  }
  stiffness_matrix_.setFromTriplets(update.begin(), update.end());
  stiffness_matrix_.makeCompressed();
}

void mpm::ProjectionSolver::apply_water_velocity_boundary_conditions_to_system() {
    if((water_velocity_constraints_.size() != stiffness_matrix_.rows()) || (water_velocity_constraints_.size() != stiffness_matrix_.cols())) {
        std::cerr << "Error file: projection_solver.ipp Line: 404" << "\n";
        std::cerr << "\t sizes of stifness matrix and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }
    if (water_velocity_constraints_.size() != force_vector_.size()) {
        std::cerr << "Error file: projection_solver.ipp Line: 409" << "\n";
        std::cerr << "\t sizes of force vector and velocity degrees of freedom are mismatched" << "\n";
        abort();
    }

 typedef Eigen::Triplet<double> T;
  std::vector<T> update;
  update.reserve(stiffness_matrix_.nonZeros() +
                 velocity_constraints_.size() * (2 * velocity_dof_ + 1));

  std::map<std::pair<unsigned, unsigned>, bool> update_vector;

  unsigned n = 0;
  for (const auto &nvelocity_bc : water_velocity_constraints_) {
    if (std::get<0>(nvelocity_bc)) {
      const auto velocity_bc = std::get<1>(nvelocity_bc);
      force_vector_ -=  velocity_bc * stiffness_matrix_.col(n);
      force_vector_[n] = velocity_bc;
      for (unsigned i = 0; i < velocity_dof_; i++) {
        // stiffness_matrix_.coeffRef(n, i) = 0;
        // stiffness_matrix_.coeffRef(i, n) = 0;
        update.push_back(T(n, i, 0));
        update.push_back(T(i, n, 0));
        update_vector.emplace(std::make_pair(std::make_pair(n, i), false));
        update_vector.emplace(std::make_pair(std::make_pair(i, n), false));
      }
      update.push_back(T(n, n, 1));
      update_vector.emplace(std::make_pair(std::make_pair(n, n), false));
    }
    n++;
  }
  // std::sort(update_vector.begin(), update_vector.end());
  for (int k = 0; k < stiffness_matrix_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness_matrix_, k);
         it; ++it) {
      if (update_vector.find(std::make_pair(it.row(), it.col())) ==
          update_vector.end())
        update.push_back(T(it.row(), it.col(), it.value()));
    }
  }
  stiffness_matrix_.setFromTriplets(update.begin(), update.end());
  stiffness_matrix_.makeCompressed();

}

void mpm::ProjectionSolver::check_for_zeroes_in_linear_system() {
  // Iterating through non-zero values
  for (int k = 0; k < stiffness_matrix_.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(stiffness_matrix_, k);
         it; ++it) {
      if (std::fabs(it.value()) < 1.E-08)
        stiffness_matrix_.coeffRef(it.row(), it.col()) = 0.;
    }
  }

    for (unsigned i = 0; i < force_vector_.size(); i++) {
        if (std::fabs(force_vector_(i)) < 1.E-08)
            force_vector_(i) = 0.;
    }
    stiffness_matrix_.makeCompressed();
}

void mpm::ProjectionSolver::assign_final_pressure_to_nodes() {
  unsigned node = 0;
  for (const auto &node_ptr : node_ptrs_) {
    double npressure = pressure_(node);
    node_ptr->update_nodal_pressure_increment(npressure);
    node_ptr->update_multimaterial_nodal_pressure_increments(0,npressure);
    node++;
  }
}

void mpm::ProjectionSolver::assign_final_acceleration_force_to_nodes() {
  Eigen::VectorXd nforce = K_S2_ * pressure_;
  unsigned node = 0;
  for (const auto node_ptr : node_ptrs_) {
    Eigen::Matrix<double,1,dim> node_force;
    if (dim == 2) {
      node_force(0) = nforce(node);
      node_force(1) = nforce(node + num_nodes_);
    }
    node_ptr->assign_KS2_force(node_force);
    node_ptr->assign_multimaterial_KS2_forces(0, node_force);
    node++;
  }
}


