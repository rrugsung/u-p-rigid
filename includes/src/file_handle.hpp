/*************************************************************************
                        Material Point Method
                       Author: Shyamini Kularathna
                         University of Cambridge
FILE: MPMSolver.hpp
**************************************************************************/
#ifndef MPM_FILE_HANDLE_H
#define MPM_FILE_HANDLE_H 
                                                                         
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> 

#include <boost/filesystem.hpp>

#include "Verify.hpp"
#include "PropertyParse.hpp"

#include "Mesh.hpp"
#include "MpmParticle.hpp"
#include "MaterialBase.hpp"
#include "ReadMaterial.hpp"
// #include "Projection.hpp"


namespace mpm {
    class FileHandle;
}

class mpm::FileHandle {

public:
  std::ifstream inputStream;
  std::ifstream nodeStream;
  std::ifstream elementStream;
  std::ifstream meshDataStream;
  std::ifstream particleStream;
  std::ifstream materialStream;
  std::ifstream iStressStream;
  std::ifstream genConstraintStream;
  std::ifstream fricConstraintStream;
  std::ifstream tractionStream;

  std::string ResultsDir;

public:
  // default constructor
  //! param[in] p path to the input files
  FileHandle(boost::filesystem::path  &p);

  // read mesh
  //! reads the input files for mesh data
  mpm::Mesh* read_mesh();

  // read particles
  //! read the input files for particle data
  mpm::MpmParticle* read_particles();

  // write output data to files
  //! param[in] step time step value
    void write_data(unsigned &step, mpm::MpmParticle* &particles, mpm::Mesh* &mesh);

    void write_para_data(unsigned &step, mpm::Mesh* &mesh);


};

#include "file_handle.ipp"

#endif
