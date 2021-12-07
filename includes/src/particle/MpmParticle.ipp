mpm::MpmParticle::MpmParticle() {
    particles_.clear();
    tp_particles_.clear();
    sp_particles_.clear();
}


void mpm::MpmParticle::read_particles(std::ifstream& p_file, std::ifstream& s_file) {
  std::string line;
  unsigned num_particles, mat_id, nphase;
  Eigen::Matrix<double, 1, dim> coordinates, spacing;
  Eigen::Matrix<double, 1, 6> init_stress;
  double volume;
  double pore_pressure;

  // read particle coordinates
  std::getline(p_file, line);
  std::istringstream in(line);
  in >> num_particles;

  // std::getline(p_file, line);
  // std::istringstream space(line);
  // for (unsigned i = 0; i < dim; i++)
  //   space >> spacing(i);

  for (unsigned i = 0; i < num_particles; i++) {
    std::getline(p_file, line);
    std::istringstream coord(line);
    for (unsigned j = 0; j < dim; j++)
      coord >> coordinates(j);
    coord >> mat_id >> volume >> nphase;
    auto particle = std::make_shared<mpm::Particle>(i, mat_id, spacing,volume,nphase);
    particle -> set_coordinates(coordinates);
    particles_.push_back(particle);
    if (nphase == 2)
        tp_particles_.push_back(particle);
    else
        sp_particles_.push_back(particle);
  }
  std::cout << "Read particles.dat" << "\n";

  // read Particle initial stress
  unsigned num_effective_stress, num_pore_pressure;
  unsigned particle_id;
  std::getline(s_file, line);
  std::istringstream is(line);
  is >> num_pore_pressure >> num_effective_stress;

  for (unsigned i = 0; i < num_pore_pressure; i++) {
      std::getline(s_file, line);
      std::istringstream pressure(line);
      pressure >> particle_id;
      pressure >> pore_pressure;
      auto particle_ptr = particles_.at(particle_id);
      particle_ptr->set_initial_pore_pressure(pore_pressure);
  }
  for (unsigned i = 0; i < num_effective_stress; i++) {
    std::getline(s_file, line);
    std::istringstream stress(line);
    stress >> particle_id;
    for (unsigned j = 0; j < 6; j++)
      stress >> init_stress(j);
    auto particle_ptr = particles_.at(particle_id);
    particle_ptr->set_initial_stress(init_stress);
  }
}


void mpm::MpmParticle::read_surface_traction(std::ifstream &traction_file) {
  unsigned num_solid_particles, num_water_particles;
  unsigned particle_id, direction;
  double surface_traction;
  std::string line;

  std::getline(traction_file, line);
  std::istringstream solid_particles(line);
  solid_particles >> num_solid_particles >> num_water_particles;

  for (unsigned i = 0; i < num_solid_particles; ++i) {
    std::getline(traction_file, line);
    std::istringstream solid_traction(line);
    solid_traction >> particle_id >> direction;
    solid_traction >> surface_traction;
    particles_.at(particle_id)->set_solid_surface_traction(direction, surface_traction);
  }

  // std::getline(traction_file, line);
  // std::istringstream water_particles(line);
  // solid_particles >> num_water_particles;
  for (unsigned i = 0; i < num_water_particles; i++) {
    std::getline(traction_file, line);
    std::istringstream water_traction(line);
    water_traction >> particle_id >> direction;
    water_traction >> surface_traction;
    particles_.at(particle_id)->set_water_surface_traction(direction, surface_traction);
  }
}

void mpm::MpmParticle::assign_material_to_particles(std::vector<mpm::material::MaterialBase*> &material_ptrs) {
    for (const auto& particle : particles_)
         particle->set_material(material_ptrs);
}

double mpm::MpmParticle::give_rigid_displ() {
    double rigid_displ = sp_particles_.at(0) -> give_rigid_displacement();
    return rigid_displ;
}

template<typename FP>
void mpm::MpmParticle::iterate_over_particles(FP function) {
  std::for_each(particles_.begin(), particles_.end(), function);
}

template<typename FP>
void mpm::MpmParticle::iterate_over_sp_particles(FP function) {
  std::for_each(sp_particles_.begin(), sp_particles_.end(), function);
}

template<typename FP>
void mpm::MpmParticle::iterate_over_tp_particles(FP function) {
  std::for_each(tp_particles_.begin(), tp_particles_.end(), function);
}


void mpm::MpmParticle::write_particle_velocity_data_to_file(std::ostream& outFile) {
    unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle Velocity Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "VECTORS Velocity float" << "\n";
    for (auto i : particles_)
        i -> write_velocity(outFile);

}




void mpm::MpmParticle::write_particle_pressure_data_to_file(std::ostream& outFile) {
    unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle Pressure Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "SCALARS Pressure float" << "\n";
    outFile << "LOOKUP_TABLE default" << "\n";
    for (auto i : particles_)
        i -> write_pressure(outFile);

}




void mpm::MpmParticle::write_particle_stress_data_to_file(std::ostream& outFile) {
    unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle Stress Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "VECTORS Stress float" << "\n";
    for (auto i : particles_)
        i -> write_stress(outFile);

}



void mpm::MpmParticle::write_particle_strain_data_to_file(std::ostream& outFile) {
    unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle Strai Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "VECTORS Strain float" << "\n";
    for (auto i : particles_)
        i -> write_strain(outFile);

}


void mpm::MpmParticle::write_particle_displacement_data_to_file(std::ostream& outFile) {
    unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle Velocity Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "VECTORS Displacement float" << "\n";
    for (auto i : particles_)
        i -> write_displacement(outFile);
}

void mpm::MpmParticle::write_particle_oed_data_to_file(std::ostream& outFile) {
    unsigned pid;
    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        pid = i->give_id();
        if(pid % 39 == 20) {
            outFile << pCord(1) << "\n";
        }
    }
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        pid = i->give_id();
        if(pid % 39 == 20) {
            i -> write_pressure(outFile);
        }
    }
}

void mpm::MpmParticle::write_particle_deviatoric_strain_data_to_file(std::ofstream& outFile) {
      unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle DeviatoricStrain Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "SCALARS DeviatoricStrain float" << "\n";
    outFile << "LOOKUP_TABLE default" << "\n";
    for (auto i : particles_)
        i -> write_deviatoric_shear_strain(outFile);
}

void mpm::MpmParticle::write_particle_plastic_strain_data_to_file(std::ofstream& outFile) {
      unsigned numOfParticles = particles_.size();

    outFile << "# vtk DataFile Version 2.0" << "\n";
    outFile << "MPM Particle PlasticStrain Data" << "\n";
    outFile << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n";
    outFile << "POINTS " << numOfParticles << " float" << "\n";

    Eigen::Matrix<double,1,dim> pCord;
    for (auto i : particles_) {
        pCord = i->give_coordinates();
        if (dim == 2)
            outFile << pCord(0) << " " << pCord(1) << " " << "0" << "\n";
        if (dim == 3)
            outFile << pCord(0) << " " << pCord(1) << " " << pCord(2) << "\n";
    }

    outFile << "CELLS " << numOfParticles << " " << 2*numOfParticles << "\n";
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1 " << i << "\n";

    outFile << "CELL_TYPES " << numOfParticles << "\n"; 
    for (unsigned i = 0; i < numOfParticles; i++)
        outFile << "1" << "\n";

    outFile << "POINT_DATA " << numOfParticles << "\n";
    outFile << "SCALARS PlasticStrain float" << "\n";
    outFile << "LOOKUP_TABLE default" << "\n";
    for (auto i : particles_)
        i -> write_equivalent_plastic_strain(outFile);
}
