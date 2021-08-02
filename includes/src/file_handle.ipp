mpm::FileHandle::FileHandle(boost::filesystem::path &p) {

    // define the input files names
    std::string inputFile = p.string() + "/inputFiles/input.dat";
    std::string nodeFile  = p.string() + "/inputFiles/node.dat" ;
    std::string elemFile  = p.string() + "/inputFiles/element.dat";
    std::string meshFile  = p.string() + "/inputFiles/meshData.dat";
    std::string particleFile = p.string() + "/inputFiles/particles.dat";
    std::string materialFile = p.string() + "/inputFiles/material.dat" ;
    std::string iStressFile = p.string() + "/inputFiles/initStress.dat";
    std::string genConstraintFile = p.string() + "/inputFiles/genCon.dat";
    std::string fricConstraintFile = p.string() + "/inputFiles/fricCon.dat";
    std::string tractionFile = p.string() + "/inputFiles/traction.dat";

    // define the input file stream for reading files
    inputStream.open(inputFile.c_str());
    nodeStream.open(nodeFile.c_str());
    elementStream.open(elemFile.c_str());
    meshDataStream.open(meshFile.c_str());
    particleStream.open(particleFile.c_str());
    materialStream.open(materialFile.c_str());
    iStressStream.open(iStressFile.c_str());
    genConstraintStream.open(genConstraintFile.c_str());
    fricConstraintStream.open(fricConstraintFile.c_str());
    tractionStream.open(tractionFile.c_str());

    // verify the file is open
    mpm::misc::VERIFY_OPEN(inputStream, inputFile);
    mpm::misc::VERIFY_OPEN(nodeStream, nodeFile);
    mpm::misc::VERIFY_OPEN(elementStream, elemFile);
    mpm::misc::VERIFY_OPEN(meshDataStream, meshFile);
    mpm::misc::VERIFY_OPEN(particleStream, particleFile);
    mpm::misc::VERIFY_OPEN(materialStream, materialFile);
    mpm::misc::VERIFY_OPEN(iStressStream, iStressFile);
    mpm::misc::VERIFY_OPEN(genConstraintStream, genConstraintFile);
    mpm::misc::VERIFY_OPEN(fricConstraintStream, fricConstraintFile);
    mpm::misc::VERIFY_OPEN(tractionStream, tractionFile);


    // create the directory for the results 
    std::string directoryOfResults = p.string() + "/Results";
    boost::filesystem::remove_all(directoryOfResults);
    boost::filesystem::path pathResults (directoryOfResults);
    if (boost::filesystem::create_directories(pathResults)) {
        std::cout << "\tResults are in " << pathResults << "\n \n";
    }
    else {
        std::cerr << "ERROR : in creating " << pathResults << "\n \n";
        abort();
    } 
    ResultsDir = pathResults.string();

    std::string line;
    while (std::getline(inputStream, line))
        mpm::misc::PARSE_INPUT_PARAMETERS(line);
    std::cout << "Input data was read" << "\n";

}


mpm::Mesh* mpm::FileHandle::read_mesh() {
  mpm::Mesh* mesh = new mpm::Mesh(meshDataStream);
  mesh->read_nodes_and_elements(nodeStream, elementStream);
  mesh->read_general_constraints(genConstraintStream);
  return mesh;

}

mpm::MpmParticle* mpm::FileHandle::read_particles() {
  mpm::MpmParticle* particle = new mpm::MpmParticle();
  particle->read_particles(particleStream, iStressStream);
  particle->read_surface_traction(tractionStream);

  // read material
  mpm::material::ReadMaterial mat(materialStream);
  std::vector<mpm::material::MaterialBase*> materials;
  materials = mat.givePtrsToMaterials();
  particle->assign_material_to_particles(materials);
  return particle;
}


void mpm::FileHandle::write_data(unsigned& step, mpm::MpmParticle* &particles, mpm::Mesh* &mesh) {
 
    std::string velocityFile = ResultsDir + "/velocity" + std::to_string(step) + ".vtk";
    std::string pressureFile = ResultsDir + "/pressure" + std::to_string(step) + ".vtk";
    std::string stressFile   = ResultsDir + "/stress" + std::to_string(step) + ".vtk";
    std::string strainFile   = ResultsDir + "/strain" + std::to_string(step) + ".vtk";
    std::string displacementFile   = ResultsDir + "/displacement" + std::to_string(step) + ".vtk";
    std::string plasticstrainFile = ResultsDir + "/plasticstrain" + std::to_string(step) + ".vtk";
    std::string deviatoricstrainFile = ResultsDir + "/deviatoricstrain" + std::to_string(step) + ".vtk";
    //std::string meshFile = ResultsDir + "/mesh" + std::to_string(step) + ".vtk";

    std::ofstream velocityOut(velocityFile.c_str()); 
    std::ofstream pressureOut(pressureFile.c_str());
    std::ofstream stressOut(stressFile.c_str());
    std::ofstream strainOut(strainFile.c_str());
    std::ofstream displacementOut(displacementFile.c_str());
    std::ofstream plasticstrainOut(plasticstrainFile.c_str());
    std::ofstream deviatoricstrainOut(deviatoricstrainFile.c_str());
    //std::ofstream meshOut(meshFile.c_str());

    particles->write_particle_velocity_data_to_file(velocityOut);
    particles->write_particle_pressure_data_to_file(pressureOut);
    particles->write_particle_stress_data_to_file(stressOut);
    particles->write_particle_strain_data_to_file(strainOut);
    particles->write_particle_displacement_data_to_file(displacementOut);
    particles->write_particle_plastic_strain_data_to_file(plasticstrainOut);
    particles->write_particle_deviatoric_strain_data_to_file(deviatoricstrainOut);
    

    velocityOut.close();
    pressureOut.close();
    stressOut.close();
    strainOut.close();
    displacementOut.close();
    plasticstrainOut.close();
    deviatoricstrainOut.close();
}

/*

//! FUNCTION: READ MPM MAIN PARAMETERS
//!          This function reads the input.dat file for main parameters
//!
void mpm::MpmSolver::ReadMainParameters() {



    gravity = mpm::misc::gravity;
    forceType = mpm::misc::forceType;
    freeSurface = mpm::misc::freeSurface;

    dt = mpm::misc::dt;
    numOfTotalSteps = mpm::misc::numOfTotalSteps;
    numOfSubSteps = mpm::misc::subSteps;

    mesh_ = new mpm::Mesh(meshDataStream);
    projectionSolver_ = new mpm::Projection(forceType);

    }


void mpm::MpmSolver::ReadModel() {

    // read mesh 
    mesh_->readNodesAndElements(nodeStream, elementStream);
    mesh_->readVelocityConstraints(velConstraintStream);
        // mesh_ -> readFrictionConstraints(fricConstraintStream);
    mesh_->computeElementCentreCoordinates();

    // read particles
    particle_.readParticles(particleStream, iStressStream);

    // read material
    mpm::material::ReadMaterial mat(materialStream);
    std::vector<mpm::material::MaterialBase*> materials;
    materials = mat.givePtrsToMaterials();

    // assign material to particles 
    particle_.assignMaterialToParticles(materials);
}


//! FUNCTION: INITIALISE ALL
//!           This function initialise the system to start the calculation
//!           at each time step
//!           This includes initialising nodes, elements, particles, matrices
//!
void mpm::MpmSolver::InitialiseAll() {

    mesh_->initialiseMesh();
    particle_.initialiseParticles();
}


//! FUNCTION: SEARCH AND LOCATE PARTICLES
//!           This function locate the particles inside the mesh and then
//!           compute the shape functions at particle locations
//!
void mpm::MpmSolver::SearchAndLocateParticles() {

    // locate particles 
    // compute shape functions at particles 
    mesh_->locateParticlesInMesh(particle_);
    particle_.computeShapeFunctionsAtParticles();

    if (freeSurface) {
        particle_.computeParticleDensityOfNodesAndElements();
        mesh_ -> defineFreeSurfaceNodes();
    }
}


//! FUNCTION: MAP DATA TO GRID
//!           This function map the particle data to the grid 
//!
void mpm::MpmSolver::MapDataToGrid() {

    particle_.mapParticleDataToGrid();

    mesh_->computeVelocityFromMomentumAtNodes();
    mesh_->computePressureFromMappedPressureAtNodes();
    mesh_->computeGravityFromMappedGravityAtNodes();

}


//! FUNCTION: SOLVE PROJECTION SYSTEM
//!           This function constructs the matrices and assemble the three main 
//!           equations in Chorin's Projection method
//!           Eq.1 : Momentum equation for solving intermediate velocity field
//!           Eq.2 : Laplace equation for solving pressure coefficient
//!           Eq.3 : Velocity Correction equation
//!
void mpm::MpmSolver::SolveProjectionSystem() {

    // build the element local matrices 
    particle_.computeElementProjectionMatrix();

    // Assembling global matrix 
    GlobalMatrixPtr globalMatrix = new GlobalMatrix(mesh_);
    globalMatrix->buildGlobalProjectionMassM();
    globalMatrix->buildGlobalProjectionViscosityM();
    globalMatrix->buildGlobalProjectionPressureGradM();
    globalMatrix->buildGlobalProjectionLaplaceM();
    globalMatrix->buildGlobalProjectionVelocityGradM();
    globalMatrix->buildGlobalProjectionProjectionCoefficientGradM();

    // Solve the system
    projectionSolver_->solveSystem(mesh_, globalMatrix, dt, particle_);

    delete globalMatrix;
}



//! FUNCTION: UPDATE THE SYSTEM
//!           This function update the particle position and velocity.
//!
void mpm::MpmSolver::UpdateTheSystem() {

    particle_.updateParticlePositionAndVelocity(dt, forceType);
    particle_.updateParticlePressure();
}

*/



