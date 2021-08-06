mpm::Node::Node(std::string& input_line, const unsigned& id) {
  nid_ = id;
  std::istringstream input(input_line);
  double third;
  input >> ncoord_(0) >> ncoord_(1);
  input >> third;

  nphase_ = 1;
  penalty_factor_ = 0;

  nsolid_mass_     = 0.;
  nwater_mass_     = 0.;
  nmixture_mass_   = 0.;
  nintrmd_mass_ = 0.;

  nsolid_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  npressure_ = 0.;
  npressure_increment_ = 0;

  nsolid_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_int_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  
  nsolid_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();  

 
  nmixture_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_int_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_force_ = Eigen::Matrix<double,1,dim>::Zero();

  KS2_force_ = Eigen::Matrix<double,1,dim>::Zero();

  nwater_damping_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_damping_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_constraints_.clear();
  nwater_constraints_.clear();
  nsolid_accel_constraints_.clear();
  nwater_accel_constraints_.clear();
  std::get<0>(npressure_constraints_) = 0;
  std::get<1>(npressure_constraints_) = 0.;
  
  undrained_status_ = false;
  moving_undrained_ = false;
  std::get<0>(undrained_node_) = 0;
  std::get<1>(undrained_node_) = 0;

  assigned_ = 0; // needs to improved
  free_pressure_ = false;
}


void mpm::Node::initialise_node() {
  nphase_ = 1;
  material_ids_.clear();

  penalty_factor_ = 0;

  nsolid_mass_     = 0.;
  nwater_mass_     = 0.;
  nmixture_mass_ = 0.;
  nintrmd_mass_ = 0.;

  nsolid_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nwater_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nmixture_masses_ = Eigen::Matrix<double,numMats,1>::Zero();
  nintrmd_masses_ = Eigen::Matrix<double,numMats,1>::Zero();

  nsolid_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_momentum_ = Eigen::Matrix<double,1,dim>::Zero();
  npressure_ = 0.;
  npressure_increment_ = 0.;

  nsolid_momenta_ = Eigen::Matrix<double,numMats,dim>::Zero();
  npressures_ = Eigen::Matrix<double,numMats,1>::Zero();
  npressure_increments_ = Eigen::Matrix<double,numMats,1>::Zero();

  nsolid_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_initial_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_final_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_initial_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nsolid_final_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();

  nsolid_int_acceleration_ = Eigen::Matrix<double,1,dim>::Zero();
  nsolid_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_velocity_ = Eigen::Matrix<double,1,dim>::Zero();

  nsolid_int_accelerations_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nsolid_int_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_int_velocities_ = Eigen::Matrix<double,numMats,dim>::Zero();
  
  nsolid_final_acceleration_ = Eigen::Matrix<double,1,dim>::Zero(); 
  nsolid_final_accelerations_ = Eigen::Matrix<double,numMats,dim>::Zero(); 
 
  nmixture_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_body_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_trac_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_int_force_ = Eigen::Matrix<double,1,dim>::Zero();
  nwater_int_force_ = Eigen::Matrix<double,1,dim>::Zero();

  nmixture_body_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_body_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nmixture_trac_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_trac_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nmixture_int_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();
  nwater_int_forces_ = Eigen::Matrix<double,numMats,dim>::Zero();

  KS2_force_ = Eigen::Matrix<double,1,dim>::Zero();

  nwater_damping_ = Eigen::Matrix<double,1,dim>::Zero();
  nmixture_damping_ = Eigen::Matrix<double,1,dim>::Zero();
    
  // if there is a moving free surface, the pressure boundary condition
  // is updated at each time step
  if(mpm::misc::freeSurface) {
    if(free_pressure_) {
      std::get<0>(npressure_constraints_) = 0;
      std::get<1>(npressure_constraints_) = 0.;
    }
  }
  if(moving_undrained_) {
    undrained_status_ = false;
    moving_undrained_ = false;
  }
  assigned_ = 0;
  free_pressure_ = false;
  pressure_status_ = true;
}

void mpm::Node::compute_nodal_initial_velocity(bool contact) {
  if (std::fabs(nsolid_mass_) > 1.0E-18)
    nsolid_initial_velocity_ = nsolid_momentum_ / nsolid_mass_;

  for (unsigned i = 0; i < dim; i++) 
    this->check_double_precision(nsolid_initial_velocity_(i));
  
  // apply velocity constraints
  if(nsolid_constraints_.size()) {
    for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
      unsigned direction = std::get<0>(nsolid_constraints_.at(j));
      nsolid_initial_velocity_(direction) = std::get<1>(nsolid_constraints_.at(j));
    }
  }

  if (contact) {
    for (unsigned i=0;  i< numMats; i++) {
      if (std::fabs(nsolid_masses_(i)) > 1.0E-18)
        nsolid_initial_velocities_.row(i) = nsolid_momenta_.row(i) / nsolid_masses_(i);

      for (unsigned j = 0; j < dim; j++)
        this->check_double_precision(nsolid_initial_velocities_(j,i));
      
      // apply velocity constraints
      if(nsolid_constraints_.size()) {
        for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
            unsigned direction = std::get<0>(nsolid_constraints_.at(j));
            nsolid_initial_velocities_(i,direction) = std::get<1>(nsolid_constraints_.at(j));
        }
      }
    }
  }
}

void mpm::Node::matrix_test() {
  if (nid_ == 399 || nid_ == 400 || nid_ == 420 || nid_ == 421) {
    std::cout << "node_id: " << nid_ << " \n";
    //std::cout << "nsolid_masses_: \n" << nsolid_masses_ << " \n";
    //std::cout << "nsolid_mass_: \n" << nsolid_mass_ << " \n";
    //std::cout << "nwater_masses_: \n" << nwater_masses_ << " \n";
    //std::cout << "nwater_mass_: \n" << nwater_mass_ << "\n";
    //std::cout << "nintrmd_masses_: \n" << nintrmd_masses_ << "\n";
    //std::cout << "mixture_masses_: \n" << nmixture_masses_ << " \n";
    //std::cout << "mixture_mass_: \n" << nmixture_mass_ << " \n\n";
    //std::cout << "nsolid_momentum: \n" << nsolid_momentum_ << "\n";
    //std::cout << "nsolid_momenta: \n" << nsolid_momenta_ << "\n";
    //std::cout << "nmixture_body_forces_: \n" << nmixture_body_forces_ << " \n";
    //std::cout << "nmixture_body_force_: \n" << nmixture_body_force_ << " \n";
    //std::cout << "nwater_body_forces_: \n" << nwater_body_forces_ << " \n";
    //std::cout << "nwater_body_force_: \n" << nwater_body_force_ << " \n";
    //std::cout << "nmixture_int_forces_: \n" << nmixture_int_forces_ << " \n";
    //std::cout << "nmixture_int_force_: \n" << nmixture_int_force_ << " \n";
    //std::cout << "nwater_int_forces_: \n" << nwater_int_forces_ << " \n";
    //std::cout << "nwater_int_force_: \n" << nwater_int_force_ << " \n\n";
    //std::cout << "nsolid_initial_velocities_: \n" << nsolid_initial_velocities_ << " \n";
    //std::cout << "nsolid_intitial_velocity_: \n" << nsolid_initial_velocity_ << " \n";  
    //std::cout << "nsolid_int_accelerations_: \n" << nsolid_int_accelerations_ << " \n";
    //std::cout << "nsolid_int_acceleration_: \n" << nsolid_int_acceleration_ << " \n";  
    //std::cout << "nsolid_int_velocities_: \n" << nsolid_int_velocities_ << " \n";
    //std::cout << "nsolid_int_velocity_: \n" << nsolid_int_velocity_ << " \n";    
    //std::cout << "nwater_int_velocities_: \n" << nwater_int_velocities_ << " \n";
    //std::cout << "nwater_int_velocity_: \n" << nwater_int_velocity_ << " \n";
    //std::cout << "nsolid_final_accelerations_: \n" << nsolid_final_accelerations_ << " \n";
    //std::cout << "nsolid_final_acceleration_: \n" << nsolid_final_acceleration_ << " \n"; 
    //std::cout << "unit_normal_vector: \n" << normal_unit_vectors_ <<"\n";
    //std::cout << "penalty_factor: \n" << penalty_factor_ << "\n";
    //std::cout << "nsolid_final_velocities_: \n" << nsolid_final_velocities_ << " \n";
    //std::cout << "nsolid_final_velocity_: \n" << nsolid_final_velocity_ << " \n\n"; 
    //std::cout << "relative_velocities: \n" << nsolid_relative_velocities_ << "\n\n";
    std::cout << "ncoord_: \n" << ncoord_ << "\n\n";
    std::cout << "=====================================================" << "\n\n";
  }
}

void mpm::Node::compute_nodal_pressure() {   
  for (unsigned i=0; i<numMats; i++){
    if (std::fabs(nwater_masses_(i)) > 1.E-18)
      npressures_(i) = npressures_(i) / nwater_masses_(i);
    this->check_double_precision(npressures_(i));
    this->apply_pressure_constraint();
  }
}

void mpm::Node::compute_intermediate_solid_acceleration(const unsigned index, const double &dt) {
  if (index == 0){
    // Compute CoM int_accelerations and velocities with the moment balance equation.
    if (nphase_ == 2){
      auto force = nmixture_body_force_ + nmixture_int_force_ + nmixture_trac_force_;
      if(std::fabs(nmixture_mass_) > 1.E-18)
        nsolid_int_acceleration_ = force / nmixture_mass_;
      nsolid_int_velocity_ = nsolid_initial_velocity_ + (dt * nsolid_int_acceleration_);  
    }

    // Update int velocity of the rigid body with the rigid_int acceleration.
    if (material_ids_.find(1) != material_ids_.end()) {
      nsolid_int_velocities_.row(1) = nsolid_initial_velocities_.row(1) + dt * nsolid_int_accelerations_.row(1);
    }
    
    // at the contact nodes, apply contact algorithm and update CoM int velocities.
    if (material_ids_.size() == 2) {
      Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_int_velocities_.row(1)-nsolid_int_velocity_;
      double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(0));
      if (std::fabs(velocity_normal) < 1E-15)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(0)*penalty_factor_;
      nsolid_int_velocity_ = nsolid_int_velocity_ + normal_correction;
    }

    // apply velocity constraints
    if(nsolid_constraints_.size()) {
        for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
            unsigned direction = std::get<0>(nsolid_constraints_.at(j));
            nsolid_int_velocity_(direction) = std::get<1>(nsolid_constraints_.at(j));
            nsolid_int_acceleration_(direction) = 0;
        }
    }
  }
  else {
    // Compute multimaterial int velocity of two-phase body using the moment balance equation.
    unsigned mat_id = index - 1;
    Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
    force = nmixture_body_forces_.row(mat_id)+nmixture_int_forces_.row(mat_id)+nmixture_trac_forces_.row(mat_id);
    if(std::fabs(nmixture_masses_(mat_id)) > 1.E-18)
      nsolid_int_accelerations_.row(mat_id) = force / nmixture_masses_(mat_id);
    nsolid_int_velocities_.row(mat_id) = nsolid_initial_velocities_.row(mat_id) + (dt * nsolid_int_accelerations_.row(mat_id));
    
    // apply velocity constraints
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_int_velocities_(mat_id,direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_int_accelerations_(mat_id,direction) = 0;
      }
    }
  }
}

void mpm::Node::compute_final_solid_acceleration(const unsigned index, const double &dt) {
  // Compute CoM velocities and accelerations at the two-phase nodes with the momentum balance equation.
  if (nphase_ == 2) {
    auto force = ((nmixture_mass_*nsolid_int_acceleration_ - KS2_force_));
    if(std::fabs(nmixture_mass_ > 1.E-18))
      nsolid_final_acceleration_ = force / nmixture_mass_;
    nsolid_final_velocity_ = nsolid_initial_velocity_ + (dt * (nsolid_final_acceleration_));
  }

  // Update final velocity of the rigid body with the rigid acceleration already assigned to the nodes.
  if (material_ids_.find(1) != material_ids_.end()) {
    nsolid_final_velocities_.row(1) = nsolid_initial_velocities_.row(1) + (dt * nsolid_final_accelerations_.row(1));
  }

  // At the contact nodes, apply contact algorithm and update the CoM velocities in normal direction.
  if (material_ids_.size() == 2) {
    Eigen::Matrix<double, 1, dim> relative_velocity = nsolid_final_velocities_.row(1) - nsolid_final_velocity_;
      double velocity_normal = relative_velocity.dot(normal_unit_vectors_.row(0));
      if (std::fabs(velocity_normal) < 1E-15)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> normal_correction = velocity_normal*normal_unit_vectors_.row(0)*penalty_factor_;
      nsolid_final_velocity_ = nsolid_final_velocity_ + normal_correction;
    //nsolid_final_velocity_ = nsolid_int_velocity_ + (dt * (nsolid_final_acceleration_ - nsolid_int_acceleration_));
  }
  // Apply velcoity constraints
  if(nsolid_constraints_.size()) {
    for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
      unsigned direction = std::get<0>(nsolid_constraints_.at(j));
      nsolid_final_velocity_(direction) = std::get<1>(nsolid_constraints_.at(j));
      nsolid_final_acceleration_(direction) = 0;
    }
  }
  
  //Compute multimaterial velocities and acclerations of the two-phase material with the momentum balance equation.
  unsigned mat_id = index - 1;
  auto force = nmixture_body_forces_.row(mat_id) + nmixture_int_forces_.row(mat_id) + nmixture_trac_forces_.row(mat_id) - KS2_forces_.row(mat_id);
  if(std::fabs(nmixture_masses_(mat_id))> 1.E-18)
    nsolid_final_accelerations_.row(mat_id) = (force / nmixture_masses_(mat_id));

  nsolid_final_velocities_.row(mat_id) = nsolid_int_velocities_.row(mat_id) + 
      (dt * (nsolid_final_accelerations_.row(mat_id) - nsolid_int_accelerations_.row(mat_id)));
    
  // apply velocity constraints
  for(auto mitr = 0; mitr < numMats; mitr++){
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_final_velocities_(mitr,direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_final_accelerations_(mitr,direction) = 0;
      }
    }
  }
}

void mpm::Node::compute_sp_final_solid_acceleration(const unsigned index, const double &dt) {
  if (index == 0){
    auto force = nmixture_body_force_ + nmixture_int_force_ + nmixture_trac_force_;
    if(std::fabs(nmixture_mass_) > 1.E-18)
      nsolid_final_acceleration_ = force / nmixture_mass_;
    nsolid_final_velocity_ = nsolid_initial_velocity_ + (dt * nsolid_final_acceleration_);
    
    // apply velocity constraints
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_final_velocity_(direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_final_acceleration_(direction) = 0;
      }
    }
  }
  else {
    unsigned mat_id = index - 1;
    Eigen::Matrix<double, 1, dim> force = Eigen::Matrix<double, 1, dim>::Zero();
    force = nmixture_body_forces_.row(mat_id)+nmixture_int_forces_.row(mat_id)+nmixture_trac_forces_.row(mat_id);
    if(std::fabs(nmixture_masses_(mat_id)) > 1.E-18)
      nsolid_final_accelerations_.row(mat_id) = force / nmixture_masses_(mat_id);
    nsolid_final_velocities_.row(mat_id) = nsolid_initial_velocities_.row(mat_id) + (dt * nsolid_final_accelerations_.row(mat_id));
    
    // apply velocity constraints
    if(nsolid_constraints_.size()) {
      for (unsigned j = 0; j < nsolid_constraints_.size(); j++) {
        unsigned direction = std::get<0>(nsolid_constraints_.at(j));
        nsolid_final_velocities_(mat_id,direction) = std::get<1>(nsolid_constraints_.at(j));
        nsolid_final_accelerations_(mat_id,direction) = 0;
      }
    }
  }
}

void mpm::Node::compute_intermediate_water_velocity(const unsigned index) {
  if (index == 0) {
    if(std::fabs(nwater_mass_) > 1.E-18) {
      nwater_int_velocity_ = (nwater_body_force_ + nwater_int_force_ + nwater_trac_force_ - (nintrmd_mass_ * nsolid_int_acceleration_)) / nwater_mass_;
      if(nwater_constraints_.size()) {
        for (unsigned i = 0; i < nwater_constraints_.size(); i++) {
          unsigned direction = std::get<0>(nwater_constraints_.at(i));
          nwater_int_velocity_(direction) = std::get<1>(nwater_constraints_.at(i));
        }
      }
      if(undrained_status_) {
        unsigned direction = std::get<0>(undrained_node_);
        nwater_int_velocity_(direction) = 0.;
      }
    }
  }
  else {
    unsigned mat_id = index - 1;
    if (std::fabs(nwater_masses_(mat_id)) > 1.E-18) {
      nwater_int_velocities_.row(mat_id) = (nwater_body_forces_.row(mat_id)+nwater_int_forces_.row(mat_id)+nwater_trac_forces_.row(mat_id) 
      - nintrmd_masses_.row(mat_id) * nsolid_int_accelerations_.row(mat_id)) / nwater_masses_(mat_id);
      if(material_ids_.size() == 2)
        nwater_int_velocities_(mat_id,1) = 0.;
      if(nwater_constraints_.size()) {
        for (unsigned i = 0; i < nwater_constraints_.size(); i++) {
          unsigned direction = std::get<0>(nwater_constraints_.at(i));
          nwater_int_velocities_(mat_id,direction) = std::get<1>(nwater_constraints_.at(i));
        }
      }
      if(undrained_status_) {
        unsigned direction = std::get<0>(undrained_node_);
        nwater_int_velocities_(mat_id,direction) = 0.;
      }
    }
  }
}

void mpm::Node::compute_multimaterial_relative_velocities(){
  for (auto mitr = 0; mitr < numMats; mitr++){
    nsolid_relative_velocities_.row(mitr) = nsolid_final_velocities_.row(mitr) - nsolid_final_velocity_;
  }
}

void mpm::Node::compute_multimaterial_normal_unit_vectors() {
  Eigen::Matrix<double,1,dim> largest_domain_gradient = Eigen::Matrix<double,1,dim>::Zero();
  double max_magnitude = 0;
  unsigned mat_id_largest = 0;
  for (unsigned i = 0; i < numMats; i++){
    Eigen::Matrix<double,1,dim> current_domain_gradient = domain_gradients_.row(i);
    if (current_domain_gradient.norm() >= max_magnitude){
      max_magnitude = current_domain_gradient.norm();
      largest_domain_gradient = current_domain_gradient;
      mat_id_largest = i;
    }
  }
  for (unsigned i = 0; i < numMats; i++) {
    Eigen::Matrix<double,1,dim> normal_unit_vector = Eigen::Matrix<double,1,dim>::Zero();
    if (largest_domain_gradient.norm() > 1.E-18)
      normal_unit_vector = largest_domain_gradient.normalized();
    if (mat_id_largest != i)
      normal_unit_vector = -1 * normal_unit_vector;
    if (nid_ == 378 || nid_ == 399 || nid_ == 420) {
      if (i == 0) {
        normal_unit_vector(0) = 0;
        normal_unit_vector(1) = 1;
      }
      else if (i == 1) {
        normal_unit_vector(0) = 0;
        normal_unit_vector(1) = -1;
      }
    }
    normal_unit_vectors_.row(i) = normal_unit_vector;
  }
}

void mpm::Node::apply_contact_mechanics(const double& dt){
  if (material_ids_.size() > 1) {
    double tolerance = 1.E-12;
    //for (auto i=0; i<numMats; i++){
      unsigned i = 0;
      Eigen::Matrix<double,1,dim> normal_unit_vector = normal_unit_vectors_.row(i);
      Eigen::Matrix<double,1,dim> relative_velocity = nsolid_relative_velocities_.row(i);
      double velocity_normal = relative_velocity.dot(normal_unit_vector);
      if (std::fabs(velocity_normal) < tolerance)
        velocity_normal = 0;
      Eigen::Matrix<double,1,dim> corrected_velocity = nsolid_final_velocities_.row(i);
      if (velocity_normal > 0){
        Eigen::Matrix<double,1,dim> corrections = Eigen::Matrix<double,1,dim>::Zero();
        double cross_product =
                relative_velocity(0) * normal_unit_vector(1) -
                relative_velocity(1) * normal_unit_vector(0);
        double mu = 0;

       // Compute the normal and tangential corrections
        Eigen::Matrix<double,1,dim> normal_correction = -velocity_normal * normal_unit_vector;
        Eigen::Matrix<double,1,dim> tangent_correction = Eigen::Matrix<double,1,dim>::Zero();
        //tangent_correction(0) = normal_unit_vector(1) * cross_product;
        //tangent_correction(1) = -normal_unit_vector(0) * cross_product;
        //tangent_correction = -mu * velocity_normal / std::abs(cross_product) *
        //                       tangent_correction;

      // Update the velocity with the computed corrections
        corrections = (normal_correction + tangent_correction);
        corrected_velocity = corrected_velocity + corrections;

        for (auto j=0; j<dim; j++){
          if (std::abs(corrected_velocity(j)) < tolerance)
              corrected_velocity(j) = 0.0;
        }
        nsolid_final_velocities_.row(i) = corrected_velocity;
        Eigen::Matrix<double,1,dim> corrected_acc = corrections / dt;
        nsolid_final_accelerations_.row(i) = nsolid_final_accelerations_.row(i) + corrected_acc;
      }
    //}
  }
}

void mpm::Node::update_mesh_configuration (double rigid_displacement, double soil_depth) {
  if (ncoord_(1) >= soil_depth)
    ncoord_(1) = ncoord_(1)+rigid_displacement;
  else
    ncoord_(1) = ncoord_(1)*((soil_depth+rigid_displacement)/soil_depth);
}

void mpm::Node::check_double_precision(double& value) {
  if (std::fabs(value) < 1.E-15)
    value = 0.;
}

void mpm::Node::print_nodal_phase() {
  std::cout << "id: " << nid_ << ", nphase: " << nphase_ << "\n";
}

void mpm::Node::apply_solid_velocity_constraint(Eigen::Matrix<double,1,dim> velocity) {
    if(nsolid_constraints_.size()) {
        for (unsigned i = 0; i < nsolid_constraints_.size(); i++) {
            unsigned direction = std::get<0>(nsolid_constraints_.at(i));
            velocity(direction) = std::get<1>(nsolid_constraints_.at(i));
        }
    }
}

void mpm::Node::apply_water_velocity_constraint(Eigen::Matrix<double,1,dim> velocity) {
  if(nwater_constraints_.size()) {
    for (unsigned i = 0; i < nwater_constraints_.size(); i++) {
      unsigned direction = std::get<0>(nwater_constraints_.at(i));
      velocity(direction) = std::get<1>(nwater_constraints_.at(i));
    }
  }
}

void mpm::Node::apply_solid_acceleration_constraint(Eigen::Matrix<double,1,dim> acceleration) {
  if(nsolid_accel_constraints_.size()) {
    for (unsigned i = 0; i < nsolid_accel_constraints_.size(); i++) {
      unsigned direction = std::get<0>(nsolid_accel_constraints_.at(i));
      acceleration(direction) = 0.0;
    }
  }
}

void mpm::Node::apply_water_acceleration_constraint(Eigen::Matrix<double,1,dim> &acceleration) {
    if(nwater_accel_constraints_.size()) {
        for (unsigned i = 0; i < nwater_accel_constraints_.size(); i++) {
            unsigned direction = std::get<0>(nwater_accel_constraints_.at(i));
            acceleration(direction) = 0.0;
        }
    }
}


void mpm::Node::apply_pressure_constraint() {
    if (std::get<0>(npressure_constraints_)) {
      npressure_ = std::get<1>(npressure_constraints_);
      for (unsigned i=0; i<numMats; i++) {
        npressures_(i) = std::get<1>(npressure_constraints_);
      }  
    }
}

void mpm::Node::apply_undrained_conditions(Eigen::Matrix<double,1,dim> water_relative_velocity) {
  if(undrained_status_) {
    unsigned direction = std::get<0>(undrained_node_);
    water_relative_velocity(direction) = 0.;
  }
}

