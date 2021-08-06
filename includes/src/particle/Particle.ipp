mpm::Particle::Particle(const unsigned& id, const unsigned& material_id, const Eigen::Matrix<double, 1, dim>& spacing,const double& vol, const unsigned& nphase) {
  id_ = id;
  mat_id_ = material_id;
  spacing_ = spacing;
  volume_ = vol;
  nphase_ = nphase;
  spacing_(0) = spacing_(1)=sqrt(volume_);
  gravity_ = Eigen::Matrix<double, 1, dim>::Zero();
  if (mpm::misc::gravity)
    gravity_(dim -1) = -9.81;

  solid_traction_.clear();
  water_traction_.clear();

  solid_mass_ = 0.;
  water_mass_ = 0.;

  displacement_ = Eigen::Matrix<double, 1, dim>::Zero();

  solid_velocity_ = Eigen::Matrix<double, 1, dim>::Zero();
  water_velocity_ = Eigen::Matrix<double, 1, dim>::Zero();

  solid_acceleration_ = Eigen::Matrix<double, 1, dim>::Zero();
  water_acceleration_ = Eigen::Matrix<double, 1, dim>::Zero();

  pore_pressure_ = 0.;
  solid_stress_ = Eigen::Matrix<double, 1, 6>::Zero();
  solid_strain_ = Eigen::Matrix<double, 1, dof>::Zero();
  solid_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
  solid_centre_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();

  element_ = NULL;

  if (dim == 2) {
    m(0) = 1;
    m(1) = 1;
    m(2) = 0;
  }
  plastic_strain_ = 0.;
}

void mpm::Particle::print_id(){
  std::cout << id_ << "\n";
}
void mpm::Particle::initialise_particle() {
    solid_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
    solid_centre_strain_rate_ = Eigen::Matrix<double, 1, dof>::Zero();
    du_ = Eigen::Matrix<double,1,dim>::Zero();
}

void mpm::Particle::set_initial_pore_pressure(double &initial_pressure) {
    pore_pressure_ = initial_pressure;
}


void mpm::Particle::set_initial_stress(const Eigen::Matrix<double, 1, 6>& initial_stress) { 
    solid_stress_ = initial_stress;
}


void mpm::Particle::set_solid_surface_traction(const unsigned &direction, const double &solidtraction) {
  std::tuple<unsigned, double> s_force = std::make_tuple(direction,solidtraction);
  solid_traction_.push_back(s_force);
}

void mpm::Particle::set_water_surface_traction(const unsigned &direction, const double& watertraction) {
  std::tuple<unsigned, double> w_force = std::make_tuple(direction,watertraction);
  water_traction_.push_back(w_force);
}


void mpm::Particle::set_material(const std::vector<mpm::material::MaterialBase*> materials) {
  material_ptrs_ = materials;
  material_ = materials.at(mat_id_);
  solid_grain_density_ = material_ -> give_density();
  water_grain_density_ = 1000.0;
  porosity_ = material_->give_porosity();
  permeability_ = material_->give_permeability();

  compute_mass();
}


void mpm::Particle::compute_mass() {
  solid_mass_ = (1 - porosity_) * solid_grain_density_ * volume_;
  water_mass_ = porosity_ * water_grain_density_ * volume_;
  if (water_mass_ < 1.0E-15)
    water_mass_ = 0.0;
}

void mpm::Particle::update_nodal_phase() {
  for (unsigned i = 0; i <numNodes; i++){
    unsigned nodal_phase = 2;
    nodes_(i)->assign_nodal_phase(nodal_phase);
  }
}

void mpm::Particle::append_material_id_to_nodes() {
  for (unsigned i = 0; i < numNodes; i++)
    nodes_(i)->append_material_id(mat_id_);
}

void mpm::Particle::map_mass_to_nodes(bool contact) {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_water_mass = water_mass_ * shape_fun_(i);
    node_mixture_mass = (solid_mass_ + water_mass_) * shape_fun_(i);
    node_intrmd_mass = water_mass_ * shape_fun_(i) * permeability_ / (porosity_ * 9.81);
    nodes_(i)->assign_nodal_masses(node_solid_mass, node_water_mass, node_mixture_mass);
    nodes_(i)->assign_nodal_intrmd_mass(node_intrmd_mass);
    if (contact) {
      nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
      nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
    }
  }
}

void mpm::Particle::map_multimaterial_domain_gradients(){
  for (unsigned i=0; i < numNodes; i++){
    Eigen::Matrix<double, 1, dim> gradient;
    gradient = volume_ * grad_shape_fun_.col(i);
    nodes_(i)->assign_multimaterial_domain_gradients(mat_id_, gradient);
  }
}

void mpm::Particle::compute_penalty_factor() {
  for (unsigned i=0; i< numNodes; i++) {
    double temp_factor;
    double factor = nodes_(i)->give_node_penalty_factor();
    Eigen::Matrix<double,1,dim> d = element_->give_element_length();
    std::set<unsigned> material_ids = nodes_(i)->material_ids_;
    if (material_ids.size() == 2) {
      Eigen::Matrix<double,1,dim> ncoord = nodes_(i)->give_node_coordinates();
      double s1 = ncoord(0) - coord_(0);
      double s2 = ncoord(1) - coord_(1);
      double s = sqrt(pow(s1,2) + pow(s2,2))-0.03536;
      if (s < 0)
        s = 0;
      //temp_factor = 1 - pow((s/d(0)),3);
      temp_factor = pow((d(0)-s)/d(0),3);
      if(temp_factor < 0) {
        temp_factor = 0;
      }
      if(temp_factor > factor){
        nodes_(i)->assign_node_penalty_factor(temp_factor);
      }
    } 
  }
}

void mpm::Particle::map_multimaterial_masses_to_nodes() {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_water_mass = water_mass_ * shape_fun_(i);
    node_mixture_mass = (solid_mass_ + water_mass_) * shape_fun_(i);
    node_intrmd_mass = water_mass_ * shape_fun_(i) * permeability_ / (porosity_ * 9.81);
    nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
    nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
  }
}

void mpm::Particle::map_sp_mass_to_nodes(bool contact) {
  double node_solid_mass = 0.;
  double node_water_mass = 0.;
  double node_mixture_mass = 0.;
  double node_intrmd_mass = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    node_solid_mass = solid_mass_ * shape_fun_(i);
    node_mixture_mass = solid_mass_ * shape_fun_(i);
    nodes_(i)->assign_nodal_masses(node_solid_mass, node_water_mass, node_mixture_mass);
    nodes_(i)->assign_nodal_intrmd_mass(node_intrmd_mass);
    if (contact)
    {
      nodes_(i)->assign_multimaterial_nodal_masses(mat_id_, node_solid_mass, node_water_mass, node_mixture_mass);
      nodes_(i)->assign_multimaterial_nodal_intrmd_masses(mat_id_, node_intrmd_mass);
    }
  }
}

void mpm::Particle::map_momentum_to_nodes(bool contact) {
  Eigen::Matrix<double, 1, dim> node_solid_momentum = Eigen::Matrix<double, 1, dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    node_solid_momentum = solid_mass_ * solid_velocity_ * shape_fun_(i);
    nodes_(i)->assign_nodal_momentum(node_solid_momentum);
    if (contact)
      nodes_(i)->assign_multimaterial_nodal_momenta(mat_id_, node_solid_momentum);
  }
}


void mpm::Particle::map_pore_pressure_to_nodes() {
  double node_pore_pressure = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    node_pore_pressure = water_mass_ * pore_pressure_ * shape_fun_(i);
    nodes_(i)->assign_nodal_multimaterial_pore_pressures(mat_id_,node_pore_pressure);
  }
}

void mpm::Particle::map_pore_pressure_from_nodes() {
  double temp_pore_pressure = 0.;
  for (unsigned i = 0; i < numNodes; i++) {
    double node_pore_pressure = nodes_(i)->give_node_pressure(mat_id_);
    temp_pore_pressure += shape_fun_(i) * node_pore_pressure;
    }
  if(std::fabs(temp_pore_pressure) < 1.0E-18)
    temp_pore_pressure = 0.;
  pore_pressure_ = temp_pore_pressure;
}


void mpm::Particle::assign_body_force_to_nodes(bool contact) {
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,1,dim> node_mixture_body_force = shape_fun_(i) * (solid_mass_ + water_mass_) * gravity_ ;
    Eigen::Matrix<double,1,dim> node_water_body_force = shape_fun_(i) * water_mass_ * gravity_ * permeability_ / (porosity_ * 9.81);
    nodes_(i)->assign_body_force(node_mixture_body_force, node_water_body_force);
    if (contact)
      nodes_(i)->assign_multimaterial_body_forces(mat_id_, node_mixture_body_force, node_water_body_force);
  }
}

void mpm::Particle::assign_sp_body_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> node_water_body_force = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i=0; i <numNodes; i++) {
    Eigen::Matrix<double,1,dim> node_mixture_body_force = shape_fun_(i) * solid_mass_ * gravity_ ;
    nodes_(i)->assign_body_force(node_mixture_body_force, node_water_body_force);
    if (contact)
      nodes_(i)->assign_multimaterial_body_forces(mat_id_, node_mixture_body_force, node_water_body_force);
  }
}


void mpm::Particle::assign_traction_force_to_nodes(bool contact, const double &time) {
  double pi = 22./7.;
  Eigen::Matrix<double,1,dim> ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned k = 0; k < water_traction_.size(); k++) {
      unsigned direction = std::get<0>(water_traction_.at(k));
      double traction = std::get<1>(water_traction_.at(k));
      nwater_traction(direction) -= beta_scalar_ * shape_fun_(i) * porosity_ * traction * volume_ / spacing_(direction);
    }
    nwater_traction = nwater_traction * permeability_ /(porosity_ * 9.81);
    ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned j = 0; j < solid_traction_.size(); j++) {
      unsigned direction = std::get<0>(solid_traction_.at(j));
      double total_traction =std::get<1>(solid_traction_.at(j));
      //if (time <= 0.02) {
      //  total_traction = total_traction * time/0.02;
      //}
      ntotal_traction(direction) += shape_fun_(i) * total_traction * volume_ / spacing_(direction);
    } 
    nodes_(i)->assign_traction_force(ntotal_traction, nwater_traction);
    if (contact)
      nodes_(i)->assign_multimaterial_traction_forces(mat_id_, ntotal_traction, nwater_traction);
  }
}

void mpm::Particle::assign_sp_traction_force_to_nodes(bool contact, const double &time) {
  double pi = 22./7.;
  Eigen::Matrix<double,1,dim> ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_traction = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {  
    ntotal_traction = Eigen::Matrix<double,1,dim>::Zero();
    for (unsigned j = 0; j < solid_traction_.size(); j++) {
      unsigned direction = std::get<0>(solid_traction_.at(j));
      double total_traction =std::get<1>(solid_traction_.at(j));
      //if (time <= 0.02) {
      //  total_traction = total_traction * time/0.02;
      //}
      ntotal_traction(direction) += shape_fun_(i) * total_traction * volume_ / spacing_(direction);
    } 
    nodes_(i)->assign_traction_force(ntotal_traction, nwater_traction);
    if (contact)
      nodes_(i)->assign_multimaterial_traction_forces(mat_id_, ntotal_traction, nwater_traction);
  }
}

void mpm::Particle::assign_internal_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> ntotal_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dof> total_temp_stress;
  if (dim == 2) {
    total_temp_stress(0) = solid_stress_(0) - beta_scalar_ * pore_pressure_;
    total_temp_stress(1) = solid_stress_(1) - beta_scalar_ * pore_pressure_;
    total_temp_stress(2) = solid_stress_(3);
  }

  Eigen::Matrix<double, dim, 1> int_force_t, int_force_w;
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    int_force_t = - volume_ * (Bi.transpose() * total_temp_stress.transpose());
    int_force_w = beta_scalar_ * volume_ * grad_shape_fun_.col(i) * porosity_ * pore_pressure_;

    ntotal_int_force = int_force_t.transpose();
    nwater_int_force = int_force_w * (permeability_ / (porosity_ * 9.81));
    nodes_(i)->assign_internal_force(ntotal_int_force, nwater_int_force);
    if (contact)
      nodes_(i) -> assign_multimaterial_internal_forces(mat_id_, ntotal_int_force, nwater_int_force);
  }
}

void mpm::Particle::assign_sp_internal_force_to_nodes(bool contact) {
  Eigen::Matrix<double,1,dim> ntotal_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> nwater_int_force = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dof> total_temp_stress;
  if (dim == 2) {
    total_temp_stress(0) = solid_stress_(0);
    total_temp_stress(1) = solid_stress_(1);
    total_temp_stress(2) = solid_stress_(3);
  }

  Eigen::Matrix<double, dim, 1> int_force_t;
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    int_force_t = - volume_ * (Bi.transpose() * total_temp_stress.transpose());
    ntotal_int_force = int_force_t.transpose();
    nodes_(i)->assign_internal_force(ntotal_int_force, nwater_int_force);
    if (contact)
      nodes_(i)->assign_multimaterial_internal_forces(mat_id_, ntotal_int_force, nwater_int_force);
  }
}

void mpm::Particle::compute_solid_strain_rate(bool contact) {
  Eigen::Matrix<double, dim, 1> node_solid_velocity;
  solid_strain_rate_ = Eigen::Matrix<double,1,dof>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,1,dim> velocity = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,dof,dim> Bi = B_.at(i);
    Eigen::Matrix<double,dof,dim> BBari = BBar_.at(i);
    if (contact)
      velocity = nodes_(i)->give_node_solid_multimaterial_initial_velocity(mat_id_);
    else
      velocity = nodes_(i)->give_node_solid_initial_velocity();
    node_solid_velocity = velocity.transpose();
    solid_strain_rate_ += (Bi * node_solid_velocity);
    //solid_strain_rate_ += (BBari * node_solid_velocity);
  }  
}


void mpm::Particle::compute_solid_strain_rate_at_elem_center(bool contact) {
  Eigen::Matrix<double, dim, 1> node_solid_velocity;
  solid_centre_strain_rate_ = Eigen::Matrix<double,1,dof>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,dof,dim> Bi_centre = BCentre_.at(i);
    Eigen::Matrix<double,1,dim> velocity = nodes_(i)->give_node_solid_multimaterial_final_velocity(mat_id_);
    node_solid_velocity = velocity.transpose();
    solid_centre_strain_rate_ += (Bi_centre * node_solid_velocity);
  }
}


void mpm::Particle::compute_solid_strain(const double dt) {
  solid_strain_ += (solid_strain_rate_ * dt);
}


void mpm::Particle::compute_solid_stress(const double dt) {
  Eigen::Matrix<double,1,dof> solid_dstrain = (dt * solid_strain_rate_);
  material_->compute_stress(solid_dstrain, solid_stress_, plastic_strain_);
}

void mpm::Particle::update_contact_velocity_and_position(const double& dt) {
  Eigen::Matrix<double,1,dim> temp_solid_final_acceleration = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_final_velocity = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_inc_velocity = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,1,dim> nsolid_acceleration = nodes_(i)->give_node_solid_multimaterial_final_acceleration(mat_id_);
    Eigen::Matrix<double,1,dim> nsolid_final_velocity = nodes_(i)->give_node_solid_multimaterial_final_velocity(mat_id_);
    Eigen::Matrix<double,1,dim> nsolid_inc_velocity = nodes_(i)->give_node_solid_multimaterial_increment_velocity(mat_id_);
    temp_solid_final_acceleration += shape_fun_(i) * nsolid_acceleration;
    temp_solid_final_velocity += shape_fun_(i) * nsolid_final_velocity;
    temp_solid_inc_velocity += shape_fun_(i) * nsolid_inc_velocity;
  }
  for (unsigned j = 0; j < dim; j++) {
    if (std::fabs(temp_solid_final_acceleration(j)) < 1.0E-18)
      temp_solid_final_acceleration(j) = 0.;
    if (std::fabs(temp_solid_final_velocity(j)) < 1.0E-18)
      temp_solid_final_velocity(j) = 0.;
  }

  coord_ += (dt * temp_solid_final_velocity);
  displacement_ += (dt * temp_solid_final_velocity);
  //solid_velocity_ += (dt * temp_solid_final_acceleration);
  solid_velocity_ = temp_solid_final_velocity;
  du_ = dt*temp_solid_final_velocity;
}

void mpm::Particle::update_velocity_and_position(const double& dt,
						 const double& time) {
  Eigen::Matrix<double,1,dim> temp_solid_final_acceleration = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_final_velocity = Eigen::Matrix<double,1,dim>::Zero();
  Eigen::Matrix<double,1,dim> temp_solid_inc_velocity = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++) {
    Eigen::Matrix<double,1,dim> nsolid_acceleration = nodes_(i)->give_node_solid_final_acceleration();
    Eigen::Matrix<double,1,dim> nsolid_final_velocity = nodes_(i)->give_node_solid_final_velocity();
    Eigen::Matrix<double,1,dim> nsolid_inc_velocity = nodes_(i)->give_node_solid_increment_velocity();
    temp_solid_final_acceleration += shape_fun_(i) * nsolid_acceleration;
    temp_solid_final_velocity += shape_fun_(i) * nsolid_final_velocity;
    temp_solid_inc_velocity += shape_fun_(i) * nsolid_inc_velocity;
  }

  for (unsigned j = 0; j < dim; j++) {
    if (std::fabs(temp_solid_final_acceleration(j)) < 1.0E-18)
      temp_solid_final_acceleration(j) = 0.;
    if (std::fabs(temp_solid_final_velocity(j)) < 1.0E-18)
      temp_solid_final_velocity(j) = 0.;
  }

  coord_ += (dt * temp_solid_final_velocity);
  displacement_ += (dt * temp_solid_final_velocity);
  solid_velocity_ += (dt * temp_solid_final_acceleration);
  //solid_velocity_ = temp_solid_final_velocity;
  //solid_velocity_ += temp_solid_inc_velocity;
}

void mpm::Particle::update_pressure() {
  double temp_pressure = 0.;
  if (nphase_ == 2) {
    for (unsigned i = 0; i < numNodes; i++) {
      double npressure = nodes_(i)->give_node_pressure_increment();
      temp_pressure += shape_fun_(i) * npressure;
    }
  }
  if (std::fabs(temp_pressure) < 1.0E-18)
    temp_pressure = 0.;
  pore_pressure_ = (beta_scalar_ * pore_pressure_) + temp_pressure;
}

void mpm::Particle::update_density_and_mass(const double &compressibility) {

}

void mpm::Particle::update_porosity(const double &dt) {
  double dvol_strain = dt * (solid_centre_strain_rate_(0) + solid_centre_strain_rate_(1));
  //volume_ = porosity_*volume_ + (1-porosity_)*volume_* (1 + dvol_strain);
  //if (nphase_ == 2) {
  //  porosity_ = porosity_/(1 + (1 - porosity_) * (1 + dvol_strain));
  //}
  porosity_ = 1 - ((1 - porosity_) / (1 + dvol_strain));

  volume_ = volume_ * (1 + dvol_strain);
  water_mass_ = porosity_ * water_grain_density_ * volume_;
}



void mpm::Particle::write_velocity(std::ostream& oFile) {
    if (dim == 2)
        oFile << solid_velocity_(0) << " " << solid_velocity_(1) << " " << "0" << "\n";
    if (dim == 3)
        oFile << solid_velocity_(0) << " " << solid_velocity_(1) << " " << solid_velocity_(2) << "\n";
}



void mpm::Particle::write_pressure(std::ostream& oFile) {
    oFile << pore_pressure_ << "\n";
}


void mpm::Particle::write_stress(std::ostream& oFile) {
    if (dim == 2)
      oFile << solid_stress_(0) << " " << solid_stress_(1) << " " << solid_stress_(3) << "\n";
    if (dim == 3)
        oFile << solid_stress_(0) << " " << solid_stress_(1) << " " << solid_stress_(2) << " " << solid_stress_(3) << " " << solid_stress_(4) << " " << solid_stress_(5) << "\n";
}


void mpm::Particle::write_strain(std::ostream& oFile) {
    if (dim == 2)
        oFile << solid_strain_(0) << " " << solid_strain_(1) << " " << solid_strain_(2) << "\n";
    if (dim == 3)
        oFile << solid_strain_(0) << " " << solid_strain_(1) << " " << solid_strain_(2) << " " << solid_strain_(3) << " " << solid_strain_(4) << " " << solid_strain_(5) << "\n";
}

void mpm::Particle::write_displacement(std::ostream& oFile) {
    if (dim == 2)
        oFile << displacement_(0) << " " << displacement_(1) << " " << 0. << "\n";
}

void mpm::Particle::write_deviatoric_shear_strain(std::ofstream& oFile) {
  oFile << dev_shear_strain_  << "\n";
}

void mpm::Particle::write_equivalent_plastic_strain(std::ofstream& oFile) {
    oFile << plastic_strain_  << "\n";
}


void mpm::Particle::compute_local_coordinates() {
    Eigen::Matrix<double,1,dim> elem_centre_coord = element_->give_element_centre_coord();
    Eigen::Matrix<double,1,dim> elem_length        = element_ -> give_element_length();
    for (unsigned i = 0; i < dim; i++) {
        xi_(i) = 2. * (coord_(i) - elem_centre_coord(i)) / elem_length(i);
        if (((std::fabs(xi_(i)) > 0.999999) && (std::fabs(xi_(i)) < 1.)) || (std::fabs(xi_(i)) > 1.))
            sign(xi_(i), 1.);
        else if ((std::fabs(xi_(i)) > 0.) && (std::fabs(xi_(i)) < 0.000001))
            xi_(i) = 0.;
    }
}


void mpm::Particle::compute_shape_functions() {
    if (dim == 2) {
        shape_fun_(0) = 0.25 * std::fabs((1 - xi_(0)) * (1 - xi_(1)));
        shape_fun_(1) = 0.25 * std::fabs((1 + xi_(0)) * (1 - xi_(1)));
        shape_fun_(2) = 0.25 * std::fabs((1 + xi_(0)) * (1 + xi_(1)));
        shape_fun_(3) = 0.25 * std::fabs((1 - xi_(0)) * (1 + xi_(1)));

        shape_fun_centre_(0) = 0.25;
        shape_fun_centre_(1) = 0.25;
        shape_fun_centre_(2) = 0.25;
        shape_fun_centre_(3) = 0.25;
    } else if (dim == 3) {
        shape_fun_(0) = fabs((1 - xi_(0)) * (1 - xi_(1)) * (1 - xi_(2))) / 8.;
        shape_fun_(1) = fabs((1 + xi_(0)) * (1 - xi_(1)) * (1 - xi_(2))) / 8.;
        shape_fun_(2) = fabs((1 + xi_(0)) * (1 - xi_(1)) * (1 + xi_(2))) / 8.;
        shape_fun_(3) = fabs((1 - xi_(0)) * (1 - xi_(1)) * (1 + xi_(2))) / 8.;
        shape_fun_(4) = fabs((1 - xi_(0)) * (1 + xi_(1)) * (1 - xi_(2))) / 8.;
        shape_fun_(5) = fabs((1 + xi_(0)) * (1 + xi_(1)) * (1 - xi_(2))) / 8.;
        shape_fun_(6) = fabs((1 + xi_(0)) * (1 + xi_(1)) * (1 + xi_(2))) / 8.;
        shape_fun_(7) = fabs((1 - xi_(0)) * (1 + xi_(1)) * (1 + xi_(2))) / 8.;
    }
    // std::cout << "Shape functions of particle: " << id_ << ":";
    // std::cout << "\t" << shape_fun_(0) << " " << shape_fun_(1) << " " << shape_fun_(2) << " " << shape_fun_(3) << "\n";
}


void mpm::Particle::compute_global_derivatives_shape_functions() {
    Eigen::Matrix<double,1,dim> L = element_ -> give_element_length();
    if (dim == 2) {
        grad_shape_fun_(0, 0) = -0.5 * (1 - xi_(1)) / L(0);
        grad_shape_fun_(0, 1) =  0.5 * (1 - xi_(1)) / L(0);
        grad_shape_fun_(0, 2) =  0.5 * (1 + xi_(1)) / L(0);
        grad_shape_fun_(0, 3) = -0.5 * (1 + xi_(1)) / L(0);

        grad_shape_fun_(1, 0) = -0.5 * (1 - xi_(0)) / L(1);
        grad_shape_fun_(1, 1) = -0.5 * (1 + xi_(0)) / L(1);
        grad_shape_fun_(1, 2) =  0.5 * (1 + xi_(0)) / L(1);
        grad_shape_fun_(1, 3) =  0.5 * (1 - xi_(0)) / L(1);
    }
    else if (dim == 3) { }
}


void mpm::Particle::compute_global_derivatives_shape_functions_at_centre() {
    Eigen::Matrix<double,1,dim> L = element_ -> give_element_length();
    if (dim == 2) {
        grad_shape_fun_centre_(0, 0) = -0.5 / L(0);
        grad_shape_fun_centre_(0, 1) =  0.5 / L(0);
        grad_shape_fun_centre_(0, 2) =  0.5 / L(0);
        grad_shape_fun_centre_(0, 3) = -0.5 / L(0);

        grad_shape_fun_centre_(1, 0) = -0.5 / L(1);
        grad_shape_fun_centre_(1, 1) = -0.5 / L(1);
        grad_shape_fun_centre_(1, 2) =  0.5 / L(1);
        grad_shape_fun_centre_(1, 3) =  0.5 / L(1);
    }
    else if (dim == 3) { }
}


void mpm::Particle::compute_B_matrix() {
    Eigen::Matrix<double,dof,dim> Bi_;
    for (unsigned i = 0; i < numNodes; i++) {
        if (dim == 2) {
            Bi_(0,0) = grad_shape_fun_(0,i); 
            Bi_(0,1) = 0.;
            Bi_(1,0) = 0.;              
            Bi_(1,1) = grad_shape_fun_(1,i);
            Bi_(2,0) = grad_shape_fun_(1,i); 
            Bi_(2,1) = grad_shape_fun_(0,i);
        }
        else if (dim == 3) {
            Bi_(0,0) = grad_shape_fun_(0,i); 
            Bi_(0,1) = 0.; 
            Bi_(0,2) = 0.;
            Bi_(1,0) = 0.; 
            Bi_(1,1) = grad_shape_fun_(1,i);
            Bi_(1,2) = 0.;
            Bi_(2,0) = 0.; 
            Bi_(2,1) = 0.;
            Bi_(2,2) = grad_shape_fun_(2,i);
            Bi_(3,0) = grad_shape_fun_(1,i); 
            Bi_(3,1) = grad_shape_fun_(0,i); 
            Bi_(3,2) = 0; 
            Bi_(4,0) = 0; 
            Bi_(4,1) = grad_shape_fun_(2,i); 
            Bi_(4,2) = grad_shape_fun_(1,i);
            Bi_(5,0) = grad_shape_fun_(2,i);
            Bi_(5,1) = 0.;
            Bi_(5,2) = grad_shape_fun_(0,i);
        }
        B_.at(i) = Bi_;
    }
}


void mpm::Particle::compute_B_matrix_at_centre() {
    Eigen::Matrix<double,dof,dim> BiCentre_;
    for (unsigned i = 0; i < numNodes; i++) {
        if (dim == 2) {
            BiCentre_(0,0) = grad_shape_fun_centre_(0,i); 
            BiCentre_(0,1) = 0.;
            BiCentre_(1,0) = 0.;              
            BiCentre_(1,1) = grad_shape_fun_centre_(1,i);
            BiCentre_(2,0) = grad_shape_fun_centre_(1,i); 
            BiCentre_(2,1) = grad_shape_fun_centre_(0,i);
        }
        BCentre_.at(i) = BiCentre_;
    }
}

void mpm::Particle::compute_BBar_matrix() {

    Eigen::Matrix<double,dof,dim> BBari_;
    double BBar00, BBar01, BBar10, BBar11;
    BBar00 = BBar01 = BBar10 = BBar11 = 0.;
    const double OneThird = 1.0/3.0;
    if (dim == 2) {
        for (unsigned i = 0; i < numNodes; i ++) {
            BBar00 = OneThird * (2*grad_shape_fun_(0,i)+grad_shape_fun_centre_(0,i));
            BBar01 = OneThird * (grad_shape_fun_centre_(1,i) - grad_shape_fun_(1,i));
            BBar10 = OneThird * (grad_shape_fun_centre_(0,i) - grad_shape_fun_(0,i));
            BBar11 = OneThird * (2*grad_shape_fun_(1,i)+grad_shape_fun_centre_(1,i));

            BBari_(0,0)= BBar00;
            BBari_(0,1)= BBar01;
            BBari_(1,0)= BBar10;
            BBari_(1,1)= BBar11;


            BBari_(2,0)= grad_shape_fun_(1,i);
            BBari_(2,1)= grad_shape_fun_(0,i);
            BBar_.at(i) = BBari_;
        }
    }
}

void mpm::Particle::compute_element_matrix_mp(const double& dt) {
  double M_ss_entry, M_ff_entry, K_bar_entry, L_entry, K_L4_entry;
  Eigen::Matrix<double, 1, dim> K_L2_entry, K_L3_entry, K_S2_entry, K_W2_entry;
  double mixture_density = (porosity_ * water_grain_density_) + ((1 - porosity_)* solid_grain_density_);
  double water_density = porosity_ * water_grain_density_;
  double laplace_constant = ((permeability_ / (water_density * 9.81))*((water_density / mixture_density)-porosity_)) - (dt/mixture_density);
  //if(id_ > 49999) {
  
  if (dim == 2) {
    for (unsigned i = 0; i < numNodes; i++) {
      for (unsigned j =0; j < numNodes; j++) {
        L_entry =  laplace_constant * volume_ * ((grad_shape_fun_(0,i) * grad_shape_fun_(0,j)) + (grad_shape_fun_(1,i) * grad_shape_fun_(1,j)));
	
	      K_L2_entry(0) = volume_ * shape_fun_(i) * grad_shape_fun_(0,j);
        K_L2_entry(1) = volume_ * shape_fun_(i) * grad_shape_fun_(1,j);
	
        K_L3_entry(0) = porosity_ * volume_ * grad_shape_fun_(0,i) * shape_fun_(j); 
        K_L3_entry(1) = porosity_ * volume_ * grad_shape_fun_(1,i) * shape_fun_(j);

        K_S2_entry(0) = volume_ * shape_fun_(i) * grad_shape_fun_(0,j);
        K_S2_entry(1) = volume_ * shape_fun_(i) * grad_shape_fun_(1,j);
	
        element_->build_Laplace_matrix_mp(i,j,L_entry);  
        element_->build_K_L2_matrix_mp(i,j,K_L2_entry);
        element_->build_K_L3_matrix_mp(i,j,K_L3_entry);
        element_->build_K_S2_matrix_mp(i,j,K_S2_entry);
      }
    }
  }
}


void mpm::Particle::sign(double& variable, double value) {
    if (variable < 0)
        variable = -value;
    if (variable > 0)
        variable = value;
}
