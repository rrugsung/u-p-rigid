mpm::Element::Element(const unsigned& id, const Eigen::Matrix<double,1,dim> &mesh_space) {
  elem_id_ = id;
  elem_length_ = mesh_space;

  neighbour_elements_.clear();
  elem_particle_density_ = 0.;
  elem_material_volume_ = 0.;
  elem_material_mass_ = 0.;
  free_ = 0;
  restrained_ = false;

  //element_matrix_ptr_ = std::make_shared<mpm::ElementMatrix>(elem_length_);
  M_SS_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_FF_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_SS_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_FF_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_BAR_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero(); 
  K_BAR_gp_= Eigen::Matrix<double, numNodes, numNodes>::Zero();
  L_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  L_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_L2_mp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L2_gp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L3_mp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L3_gp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L4_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_L4_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_S2_mp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_S2_gp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_W2_mp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_W2_gp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
}

void mpm::Element::initialise_element() {
  elem_particle_density_ = 0.;
  elem_material_volume_ = 0.;
  elem_material_mass_ = 0.;
  free_surface_nodes_.clear();
  free_ = 0;
  restrained_ = false;

  M_SS_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_FF_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_SS_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  M_FF_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();

  K_BAR_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero(); 
  K_BAR_gp_= Eigen::Matrix<double, numNodes, numNodes>::Zero();

  L_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  L_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();

  K_L2_mp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L2_gp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L3_mp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L3_gp_ = Eigen::Matrix<double, numNodes, dim*numNodes>::Zero();
  K_L4_mp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();
  K_L4_gp_ = Eigen::Matrix<double, numNodes, numNodes>::Zero();

  K_S2_mp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_S2_gp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_W2_mp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
  K_W2_gp_ = Eigen::Matrix<double, dim*numNodes, numNodes>::Zero();
 }


void mpm::Element::compute_centre_coordinates() {
  elem_centre_coord_ = Eigen::Matrix<double,1,dim>::Zero();
  for (unsigned i = 0; i < numNodes; i++)
    elem_centre_coord_ += elem_nodes_ptrs_(i) -> give_node_coordinates();
  elem_centre_coord_ /= numNodes;

  // algorithm for element length : Only for 2D
  Eigen::Matrix<double,1,dim> first_coord = elem_nodes_ptrs_(0)->give_node_coordinates();
  Eigen::Matrix<double,1,dim> corner_coord = elem_nodes_ptrs_(2)->give_node_coordinates();
  for (unsigned i = 0; i < dim; i++)
    elem_length_(i) = std::fabs(corner_coord(i) - first_coord(i)); 

  elem_volume_ = 0;
  if (dim == 2)
    elem_volume_ = elem_length_(0) * elem_length_(1) * 1;
  if (dim == 3)
    elem_volume_ = elem_length_(0) * elem_length_(1) * elem_length_(2);
}


void mpm::Element::find_free_surface_nodes() {
  free_ = 0;
  free_surface_nodes_.clear();
  unsigned nid = 2;
  auto nptr = this->give_element_node_ptr_at(nid);
  Eigen::Matrix<double,1,dim> ncoord = nptr->give_node_coordinates();
  //if(ncoord(0) > 0.5) {
  for (const auto &neighbour_ptr : neighbour_elements_) {
    if ((neighbour_ptr->give_ratio_particle_volume_to_element_volume()) < alpha) {
	    free_ = 1;
	    this->insert_to_free_surface_nodes(neighbour_ptr);
      }
    }
    //}
}


void mpm::Element::insert_to_free_surface_nodes(const mpm::Element* elem_ptr) {
  Eigen::Matrix<unsigned,1,numNodes> node_ids = elem_ptr->give_element_node_ids();
  for (unsigned i = 0; i < numNodes; i++) {
    unsigned element_node_id = elem_nodes_id_(i);
    if((node_ids.array() == element_node_id).any()) 
      free_surface_nodes_.insert(elem_nodes_ptrs_(i));
  }
}

void mpm::Element::build_M_matrix_mp(const unsigned& row, const unsigned& column, const double& value_ss, const double& value_ff) {
  M_SS_mp_(row, column) += (4 * value_ss / elem_material_volume_);
  M_FF_mp_(row, column) += (4 * value_ff / elem_material_volume_);
}

void mpm::Element::build_K_bar_matrix_mp(const unsigned& row, const unsigned& column, const double& value) {
  K_BAR_mp_(row, column) += (4 * value / elem_material_volume_);
}

void mpm::Element::build_Laplace_matrix_mp(const unsigned& row, const unsigned& column, const double& value) {
  L_mp_(row, column) += value;
}

void mpm::Element::build_K_L2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value) {
  K_L2_mp_(row, column) += value(0);
  K_L2_mp_(row, (numNodes + column)) +=  value(1);
}

void mpm::Element::build_K_L3_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value) {
  K_L3_mp_(row, column) += value(0);
  K_L3_mp_(row, (numNodes + column)) += value(1);
}

void mpm::Element::build_K_L4_matrix_mp(const unsigned& row, const unsigned& column, const double &value) {
  K_L4_mp_(row, column) += (4 * value / elem_material_volume_);
}

void mpm::Element::build_K_S2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value) {
  K_S2_mp_(row, column) += value(0);
  K_S2_mp_((numNodes + row), column) += value(1);
}

void mpm::Element::build_K_W2_matrix_mp(const unsigned& row, const unsigned& column, const Eigen::Matrix<double,1,dim>& value) {
  K_W2_mp_(row, column) += (4 * value(0) / elem_material_volume_);
  K_W2_mp_((numNodes + row), column) += (4 * value(1) / elem_material_volume_);
}


void mpm::Element::compute_element_matrix_gp(const double& permeability) {
  // Eigen::Matrix<double,1,numNodes> nodal_solid_porosity; 
  // Eigen::Matrix<double,1,numNodes> nodal_water_porosity;
  // unsigned node = 0;
  // for(unsigned i = 0; i < numNodes; i++) {
  //   nodal_solid_porosity(i) = elem_nodes_ptrs_(i)->give_node_solid_porosity();
  //   nodal_water_porosity(i) = elem_nodes_ptrs_(i)->give_node_water_porosity();
  // }
  // element_matrix_ptr_->assign_porosity_to_gauss_points(nodal_solid_porosity, nodal_water_porosity);
  // M_SS_gp_ = element_matrix_ptr_->M_ss_matrix_gp();
  // M_FF_gp_ = element_matrix_ptr_->M_ff_matrix_gp();
  // K_BAR_gp_ = element_matrix_ptr_->K_bar_matrix_gp(permeability);
  // L_gp_ = element_matrix_ptr_->Laplace_matrix_gp();
  // K_L2_gp_ = element_matrix_ptr_->K_L2_matrix_gp();
  // K_L3_gp_ = element_matrix_ptr_->K_L3_matrix_gp();
  // K_L4_gp_ = element_matrix_ptr_->K_L4_matrix_gp();
  // K_S2_gp_ = element_matrix_ptr_->K_S2_matrix_gp();
  // K_W2_gp_ = element_matrix_ptr_->K_W2_matrix_gp();
}
