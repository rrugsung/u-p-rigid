template<typename FP>
void mpm::Mesh::iterate_over_elements(FP function) const {
  std::for_each(elements_.begin(), elements_.end(), function);
  return;
}

template<typename FP>
void mpm::Mesh::iterate_over_nodes(FP function) const {
  std::for_each(nodes_.begin(), nodes_.end(), function);
  return;
}

template<typename FP>
void mpm::Mesh::iterate_over_elements_of_p(FP function) const {
  std::for_each(p_element_set_.begin(), p_element_set_.end(), function);
  return;
}

template<typename FP>
void mpm::Mesh::iterate_over_nodes_of_p(FP function) const {
  std::for_each(p_node_set_.begin(), p_node_set_.end(), function);
  return;
}

void mpm::Mesh::write_mesh_data_to_file(std::ostream& outFile) {
    unsigned numofnodes = nodes_.size();
    unsigned numofelements = elements_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Mesh Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numofnodes << " float " << "\n";

    Eigen::Matrix<double,1,dim> nCord;
    for (auto i : nodes_){
        nCord = i->give_node_coordinates();
        if (dim == 2)
            outFile << nCord(0) << " " << nCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << nCord(0) << " " << nCord(1) << " " << nCord(2) << "\n";
    }

    outFile << "CELLS " << numofelements << " " << 5*numofelements << "\n";
    Eigen::Matrix<unsigned,1,numNodes> nIds;
    for (auto i : elements_) {
        nIds = i->give_element_node_ids();
        outFile << numNodes << " " << nIds(0) << " " << nIds(1) 
                << " " << nIds(2) << " " << nIds(3) << "\n";
    }

    outFile << "CELL_TYPES " << numofelements << "\n";
    for (unsigned i=0; i<numofelements; i++)
        outFile << "9" << "\n";
}

void mpm::Mesh::initialise_mesh() {
  p_element_set_.clear();
  solver_element_set_.clear();
  p_node_set_.clear();
  part_fill_elems_.clear();
  free_node_set_.clear();
  for (const auto& elem : elements_)
    elem->initialise_element();
  for (const auto& node : nodes_)
    node->initialise_node();
}

void mpm::Mesh::free_memory() {
  for (auto &ePtr : elements_)
    delete ePtr;
  for (auto &nPtr : nodes_)
    delete nPtr;
  elements_.clear();
  nodes_.clear();
  p_element_set_.clear();
  solver_element_set_.clear();
  p_node_set_.clear();
  part_fill_elems_.clear();
  free_node_set_.clear();
}


mpm::Mesh::Mesh(std::ifstream& mesh_data_file) {
    elements_.clear();
    nodes_.clear();
    p_element_set_.clear();
    p_node_set_.clear();
    part_fill_elems_.clear();
    free_node_set_.clear();

    std::string line;
    std::getline(mesh_data_file, line);
    std::istringstream input_mesh_space(line);
    for (unsigned i = 0; i < dim; i++)
        input_mesh_space >> mesh_spacing_(i);

    std::getline(mesh_data_file, line);
    std::istringstream input_num_elem(line);
    for (unsigned i = 0; i < dim; i++)
        input_num_elem >> num_elements_(i);

    std::getline(mesh_data_file, line);
    std::istringstream input_corner_elem(line);
    for (unsigned i = 0; i < numCorners; i++) 
        input_corner_elem >> corner_elements_(i);

    std::getline(mesh_data_file, line);
    std::istringstream input_corner_node(line);
    for (unsigned i = 0; i < numCorners; i++)
        input_corner_node >> corner_nodes_(i);
    std::cout << "Read meshData.dat" << "\n";
}


void mpm::Mesh::read_nodes_and_elements(std::ifstream& node_file, std::ifstream& elem_file) {
    std::string line;
    // read Nodes
    unsigned num_nodes;
    if(!std::getline(node_file, line)) {
        std::cerr << "ERROR: reading file" << "\n";
        abort();
    } 
    std::istringstream input_num_nodes(line);
    input_num_nodes >> num_nodes;

    for (unsigned i = 0; i < num_nodes; i++) {
        if (!std::getline(node_file, line)) {
            std::cerr << "ERROR: reading file" << "\n";
            abort();
        }
        mpm::Node* node = new mpm::Node(line,i);
        nodes_.push_back(node);
    }
    std::cout << "Read node.dat" << "\n";

    // read Elements
    unsigned num_elements;
    if(!std::getline(elem_file, line)) {
        std::cerr << "ERROR: reading file" << "\n";
        abort();
    }
    std::istringstream input_num_elems(line);
    input_num_elems >> num_elements;

    for (unsigned i = 0; i < num_elements; i++) {
        mpm::Element* element = new mpm::Element(i, mesh_spacing_);
        if (!std::getline(elem_file, line)) {
          std::cerr << "ERROR: reading file" << "\n";
          abort();
        }
        std::istringstream input_elem_nodes(line);
        for (unsigned j = 0; j < numNodes; j++) {
            unsigned node_id;
            input_elem_nodes >> node_id;
            element->set_element_nodes(j, node_id, nodes_.at(node_id));
        }
        elements_.push_back(element);
    }

    // read the neighbour elements
    for (unsigned i = 0; i < num_elements; i++) {
        unsigned num_neighbour_elems;
        if (!std::getline(elem_file, line)) {
          std::cerr << "ERROR: reading file" << "\n";
          abort();
        }
        std::istringstream input_neighbours(line);
        input_neighbours >> num_neighbour_elems;
        unsigned elem_id;
        for (unsigned j = 0; j < num_neighbour_elems; j++) {
            input_neighbours >> elem_id;
            elements_.at(i)->set_neighbour_elements(elements_.at(elem_id));
        } 
    }
    std::cout << "Read element.dat" << "\n";
    first_node_coord_ = nodes_.at(corner_nodes_(0)) -> give_node_coordinates();
    last_node_coord_  = nodes_.at(corner_nodes_(3)) -> give_node_coordinates();
}


void mpm::Mesh::read_general_constraints(std::ifstream& vel_con_file) {
    std::string line;
    unsigned num_undrained_nodes, undrained_direction; 
    int undrained_normal;
    unsigned num_solid_con_nodes, num_water_con_nodes, num_pres_con_nodes;
    unsigned node_id, direction;
    double value;

    std::getline(vel_con_file, line);
    std::istringstream input_undrained_nums(line);
    input_undrained_nums >> num_undrained_nodes ;
    std::cout << num_undrained_nodes << "\n";
    for (unsigned i = 0; i < num_undrained_nodes; i++) {
        std::getline(vel_con_file, line);
        std::istringstream input_undrained_node(line);
        input_undrained_node >> node_id >> undrained_direction;
        input_undrained_node >> undrained_normal;
        nodes_.at(node_id) -> set_undrained_nodes(undrained_direction, undrained_normal); 
    }


    std::getline(vel_con_file, line);
    std::istringstream input_solid_nums(line);
    input_solid_nums >> num_solid_con_nodes ;
    std::cout << "solid velocity constraints: " << num_solid_con_nodes << "\n";
    for (unsigned i = 0; i < num_solid_con_nodes; i++) {
        std::getline(vel_con_file, line);
        std::istringstream input_solid_constraint(line);
        input_solid_constraint >> node_id >> direction;
        input_solid_constraint >> value;
        nodes_.at(node_id) -> set_node_solid_velocity_constraints(direction, value); 
    }

    std::getline(vel_con_file, line);
    std::istringstream input_water_nums(line);
    input_water_nums >> num_water_con_nodes;
    std::cout << "water velocity constraints: " << num_water_con_nodes << "\n";
    for (unsigned i = 0; i < num_water_con_nodes; i++) {
        std::getline(vel_con_file, line);
        std::istringstream input_water_constraint(line);
        input_water_constraint >> node_id >> direction;
        input_water_constraint >> value;
        nodes_.at(node_id) -> set_node_water_velocity_constraints(direction, value); 
    }

    std::getline(vel_con_file, line);
    std::istringstream input_pressure_nums(line);
    input_pressure_nums >> num_pres_con_nodes;
    std::cout << "Pressure boundary nodes: " << num_pres_con_nodes << "\n";  
    for (unsigned i = 0; i < num_pres_con_nodes; i++) {
        std::getline(vel_con_file, line);
        std::istringstream input_pres_con(line);
        input_pres_con >> node_id >> value;
        nodes_.at(node_id)->set_node_pressure_constraints(value);
    }
    std::cout << "Read genCon.dat" << "\n";
}


void mpm::Mesh::read_friction_constraints(std::ifstream& fric_con_file) {

}


void mpm::Mesh::locate_particles_in_mesh(mpm::MpmParticle* &particle_set) {
    unsigned num_particles = particle_set->number_of_particles();
    Eigen::Matrix<int, 1 , dim> elem_grid;
    for (unsigned i = 0; i < num_particles; i++) {
        Eigen::Matrix<double, 1 , dim> p_coords = particle_set->particle_coordinates(i);
        for (unsigned j = 0; j < dim; j++)
            elem_grid(j) = std::fabs((p_coords(j) - first_node_coord_(j)) / mesh_spacing_(j));
        check_particle_is_inside_mesh(elem_grid, i);
        unsigned elem_id;
        if (dim == 2)
            elem_id = std::fabs(elem_grid(0) + (num_elements_(0) * elem_grid(1)));
        if (elem_id > corner_elements_(3) || elem_id < corner_elements_(0)) {
            std::cerr << "ERROR: in computing element id for particle " << i << "\n";
            abort();
        }
        std::shared_ptr<mpm::Particle> pPtr = particle_set->pointer_to_particle(i);
        set_elements_and_nodes_of_particles(elem_id, pPtr);
    }
}

void mpm::Mesh::check_particle_is_inside_mesh(Eigen::Matrix<int, 1 , dim> &e_grid, unsigned &p_id) {
  for (unsigned i = 0; i < dim; i++) {
    if (e_grid(i) < 0 || e_grid (i) >= num_elements_(i)) {
      std::cerr << "ERROR: Particle " << p_id << " is out of mesh" << "\n";
      abort();
    }
  }
}


void mpm::Mesh::set_elements_and_nodes_of_particles(unsigned &element_id, std::shared_ptr<mpm::Particle> &particle_ptr) {
    mpm::Element* elem_ptr = elements_.at(element_id);
    p_element_set_.insert(elem_ptr);
    Eigen::Matrix<mpm::Node*, 1, numNodes> node_ptrs = elem_ptr->give_element_node_ptrs();
    for (unsigned i = 0; i < numNodes; i++) 
        p_node_set_.insert(node_ptrs(i));
    particle_ptr->set_element_and_nodes(elem_ptr, node_ptrs);

    double p_volume = particle_ptr->give_volume();
    double p_density = particle_ptr->give_density();
    if(particle_ptr->give_id()>74999)
       p_volume = 0.001;
    
    elem_ptr->add_to_element_particle_density(p_volume, p_density);
    if(particle_ptr->give_id() < 74999)
       solver_element_set_.insert(elem_ptr);
}


void mpm::Mesh::find_partially_filled_elements() {
    part_fill_elems_.clear();
    for (const auto &pelem_ptr : p_element_set_) {
        double e_volume = pelem_ptr->give_element_volume();
        double e_mat_volume = pelem_ptr->give_material_volume_in_element();
        if(std::fabs(e_mat_volume) < (beta * std::fabs(e_volume))){
	//if (pelem_ptr->give_ratio_particle_mass_to_element_mass() < beta){
            part_fill_elems_.push_back(pelem_ptr);
	}
    }
}

void mpm::Mesh::prescribe_pressure_at_free_surface(const double &pressure) {
    std::vector<mpm::Node*> vec_node_ptrs;
    for (const auto &p_elem_ptr : p_element_set_) {
        vec_node_ptrs.clear();
        p_elem_ptr->give_element_free_surface_nodes(vec_node_ptrs);
	// if((p_elem_ptr->give_ratio_particle_volume_to_element_volume()) <0.125) {
	//   auto enodes = p_elem_ptr->give_element_node_ptrs();
	//   for (unsigned i = 0; i < numNodes; i++)
	//     enodes(i)->set_apply_nodal_pressure_status();
	// }
        if (vec_node_ptrs.size()) {
            for (const auto &n_ptr : vec_node_ptrs)
                free_node_set_.insert(n_ptr);
        }
    }
    for (const auto &n_ptr : free_node_set_) {
        n_ptr->set_free_node_pressure_constraints(pressure);
    }
}


void mpm::Mesh::give_elements_with_particles_above_cutoff(std::vector<mpm::Element*>& vec_of_elem, const double &cut_off) {
  for (const auto &elem_ptr : p_element_set_) {
    //for (const auto &elem_ptr : solver_element_set_) {
        double mat_volume_to_elem_volume = elem_ptr->give_ratio_particle_volume_to_element_volume();
	    double mat_mass_to_elem_mass = elem_ptr->give_ratio_particle_mass_to_element_mass();
        Eigen::Matrix<mpm::Node*, 1, numNodes> temp_nodes = elem_ptr->give_element_node_ptrs();
        unsigned temp = 0;

        for (unsigned i = 0; i < numNodes; i++) {
            unsigned node_phase_ = temp_nodes(i) -> give_node_phase();
            unsigned nid = temp_nodes(i) -> give_id();
            temp = temp + node_phase_;
        }
        
        if (mat_volume_to_elem_volume > cut_off){
            if (temp == 8){
                vec_of_elem.push_back(elem_ptr);
            }
	    }
    }
}

void mpm::Mesh::compute_rigid_body_int_acceleration() {
    double rigid_mass = 0;
    Eigen::Matrix<double,1,dim> mixture_body_force = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,1,dim> mixture_int_force = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,1,dim> mixture_trac_force = Eigen::Matrix<double,1,dim>::Zero();
    for (const auto &n_ptr : p_node_set_) {
        std::set<unsigned> material_ids = n_ptr->material_ids_;
        if (material_ids.find(1) != material_ids.end()) {
            rigid_mass += n_ptr-> give_node_multimaterial_mixture_masses(1);
            mixture_body_force = mixture_body_force + n_ptr->give_node_multimaterial_mixture_body_forces(1)
                                + n_ptr->give_node_multimaterial_mixture_body_forces(0);
            mixture_trac_force = mixture_trac_force + n_ptr->give_node_multimaterial_mixture_trac_forces(1);
            Eigen::Matrix<double,1,dim> internal_force = n_ptr->give_node_multimaterial_mixture_internal_forces(0);
            Eigen::Matrix<double,1,dim> normal_vector = Eigen::Matrix<double,1,dim>::Zero();
            //normal vector is set to (0,1) as the simulation is symmetric
            normal_vector(1) = 1.;
            double force_normal = internal_force.dot(normal_vector);
            Eigen::Matrix<double,1,dim> int_normal_force = force_normal*normal_vector;
            mixture_int_force = mixture_int_force + int_normal_force;
        }
    }
    // compute rigid body acceleration
    Eigen::Matrix<double,1,dim> rigid_acceleration = Eigen::Matrix<double,1,dim>::Zero();
    if (rigid_mass > 1.E-15) {
        rigid_acceleration = (mixture_trac_force + mixture_body_force + mixture_int_force) / rigid_mass;
        //std::cout << "mixture_body_force: \n" << mixture_body_force << "\n";
        //std::cout << "mixture_int_force: \n" << mixture_int_force << "\n";
        //std::cout << "mixture_trac_force: \n" << mixture_trac_force << "\n";
        //std::cout << "int_rigid_acc_: \n" << rigid_acceleration << "\n";
    }
    for (const auto &n_ptr : p_node_set_) {
        std::set<unsigned> material_ids = n_ptr->material_ids_;
        if (material_ids.find(1) != material_ids.end()){
            n_ptr->assign_multimaterial_nsolid_int_accelerations(1,rigid_acceleration);
            //n_ptr->assign_nsolid_int_acceleration(rigid_acceleration);
        }
    }
}

void mpm::Mesh::compute_rigid_body_final_acceleration() {
    double rigid_mass = 0;
    Eigen::Matrix<double,1,dim> mixture_body_force = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,1,dim> mixture_int_force = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,1,dim> mixture_trac_force = Eigen::Matrix<double,1,dim>::Zero();
    Eigen::Matrix<double,1,dim> KS2_force = Eigen::Matrix<double,1,dim>::Zero();
    for (const auto &n_ptr : p_node_set_) {
        std::set<unsigned> material_ids = n_ptr->material_ids_;
        if (material_ids.find(1) != material_ids.end()) {
            rigid_mass += n_ptr-> give_node_multimaterial_mixture_masses(1);
            mixture_body_force = mixture_body_force + n_ptr->give_node_multimaterial_mixture_body_forces(1)
                                + n_ptr->give_node_multimaterial_mixture_body_forces(0);
            mixture_trac_force = mixture_trac_force + n_ptr->give_node_multimaterial_mixture_trac_forces(1);
            KS2_force = n_ptr->give_node_KS2_forces();
            Eigen::Matrix<double,1,dim> internal_force = n_ptr->give_node_multimaterial_mixture_internal_forces(0)
                                        - KS2_force;
            Eigen::Matrix<double,1,dim> normal_vector = Eigen::Matrix<double,1,dim>::Zero();
            //normal vector is set to (0,1) as the simulation is symmetric
            normal_vector(1) = 1.;
            double force_normal = internal_force.dot(normal_vector);
            Eigen::Matrix<double,1,dim> int_normal_force = force_normal*normal_vector;
            mixture_int_force = mixture_int_force + int_normal_force;
            
        }
    }
    // compute rigid body acceleration
    Eigen::Matrix<double,1,dim> rigid_acceleration = Eigen::Matrix<double,1,dim>::Zero();
    if (rigid_mass > 1.E-15) {
        rigid_acceleration = (mixture_trac_force + mixture_body_force + mixture_int_force) / rigid_mass;
        //std::cout << "mixture_body_force: \n" << mixture_body_force << "\n";
        //std::cout << "mixture_trac_force: \n" << mixture_trac_force << "\n";
        //std::cout << "mixture_int_force: \n" << mixture_int_force << "\n";
        //std::cout << "KS2_Force: \n" << KS2_force << "\n";
        //std::cout << "final_rigid_acc_: \n" << rigid_acceleration << "\n";
    }
    for (const auto &n_ptr : p_node_set_) {
        std::set<unsigned> material_ids = n_ptr->material_ids_;
        if (material_ids.find(1) != material_ids.end()){
            n_ptr->assign_multimaterial_final_accelerations(1,rigid_acceleration);
            //n_ptr->assign_nsolid_final_acceleration(rigid_acceleration);
        }
    }
}