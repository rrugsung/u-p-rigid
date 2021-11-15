/***************************************************************************
                          Material Point Method
                           Shyamini Kularathna
                         University Of Cambridge
File: Node.hpp
****************************************************************************/
#ifndef MPM_NODE_H
#define MPM_NODE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tuple>

#include <eigen3/Eigen/Dense>

#include "PropertyParse.hpp"
#include "Constants.hpp"
#include "ReadMaterial.hpp"

namespace mpm {
    class Node;
}

class mpm::Node {

protected:
    static const unsigned dim = mpm::constants::DIM;
    static const unsigned numNodes = mpm::constants::NUMNODES;
    static const unsigned numMats = mpm::constants::NUMBODIES;

public:
  // constructor
  //! param[in] input_line
  //! param[in] id node id
  Node(std::string& input_line, const unsigned& id);

  // initialise node
  void initialise_node();

  // set node solid velocity constraints
  //! param[in] dir direction of the velocity constraint
  //! param[in] value value of the velocity constraint
  void set_node_solid_velocity_constraints(const unsigned& dir,
					   const double value) {
    auto t = std::make_tuple(dir, value);
    nsolid_constraints_.push_back(t);
    auto at = std::make_tuple(dir,0.);
    nsolid_accel_constraints_.emplace_back(at);
    fixed_boundary_ = true;
  }

  // set node water velocity constraints
  //! param[in] dir direction of the velocity constraint
  //! param[in] value value of the velocity constraint
  void set_node_water_velocity_constraints(const unsigned& dir,
					   const double value) {
    auto t = std::make_tuple(dir, value);
    nwater_constraints_.push_back(t);
    auto at = std::make_tuple(dir,0.);
    nwater_accel_constraints_.emplace_back(at);
  }

  // set node pressure constraints
  //! param[in] value value of the prescribed pressure at the node
  void set_node_pressure_constraints(const double& value) {
    std::get<0>(npressure_constraints_) = 1;
    std::get<1>(npressure_constraints_) = value;
  }

  void set_free_node_pressure_constraints(const double& value) {
    std::get<0>(npressure_constraints_) = 1;
    std::get<1>(npressure_constraints_) = value;
    free_pressure_ = true;
  }

  // set boundary nodes and its normal direction
  //! param[in] status boundary node or not
  //! param[in] norm_direction noraml direction to the boundary
  void set_undrained_nodes(const unsigned& direction, const int& normal) {
    undrained_status_ = 1;
    std::get<0>(undrained_node_) = direction;
    std::get<1>(undrained_node_) = normal;
  }
 
  // change nodal phase
  void assign_nodal_phase(unsigned nodal_phase){
    nphase_ = nodal_phase;
  }
  
  void assign_node_penalty_factor(double factor){
    penalty_factor_ = factor;
  }

  //! Add material id from material points to list of materials in materials_
  //! \param[in] id Material id to be stored at the node
  void append_material_id(unsigned id) {
    material_ids_.emplace(id);
  };

  // print nodal phase
  void print_nodal_phase();
  
  // assign nodal lumped mass from particles
  //! param[in] mass lumped mass from particles : OK
  void assign_nodal_masses(const double& solid_mass, const double& water_mass,
			   const double& mixture_mass) {
    nsolid_mass_ += solid_mass;
    nwater_mass_ += water_mass;
    nmixture_mass_ += mixture_mass;
  }

  // assig nodal multimaterial domain gradient
  //! param[in] material id, domain gradient
  void assign_multimaterial_domain_gradients(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& gradient){
    domain_gradients_.row(mat_id) += gradient;
  }

  // compute multimaterial unit normal vector
  void compute_multimaterial_normal_unit_vectors();

  // assign multimaterial_nodal lumped mass from particles
  //! param[in] mass lumped mass from particles : OK
  void assign_multimaterial_nodal_masses(const unsigned mat_id, const double& solid_mass, const double& water_mass,
			   const double& mixture_mass) {
    nsolid_masses_(mat_id) += solid_mass;
    nwater_masses_(mat_id) += water_mass;
    nmixture_masses_(mat_id) += mixture_mass;
  }

  // assign multimaterial solid intermediate accelerations
  //! param[in] mat_id, nsolid_int_accelerations_
  void assign_multimaterial_nsolid_initial_velocities(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& velocity){
    nsolid_initial_velocities_.row(mat_id) = velocity;
  }

  // assign multimaterial solid intermediate accelerations
  //! param[in] mat_id, nsolid_int_accelerations_
  void assign_multimaterial_nsolid_int_accelerations(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& acceleration){
    nsolid_int_accelerations_.row(mat_id) = acceleration;
  }

  // assign multimaterial solid final accelerations
  //! param[in] mat_id, nsolid_final_accelerations_
  void assign_multimaterial_final_accelerations(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& acceleration){
    nsolid_final_accelerations_.row(mat_id) = acceleration;
  }

  // assign solid final acceleration
  //! param[in] nsolid_final_acceleration_
  void assign_nsolid_final_acceleration(const Eigen::Matrix<double,1,dim>& acceleration){
    nsolid_final_acceleration_ = acceleration;
  }

  // assign solid intermediate acceleration
  //! param[in] nsolid_int_acceleration_
  void assign_nsolid_int_acceleration(const Eigen::Matrix<double,1,dim>& acceleration){
    nsolid_int_acceleration_ = acceleration;
  }

  // assign nodal momentum from particles
  //! param[in] momentum momentum from particles : OK
  void assign_nodal_momentum(const Eigen::Matrix<double,1,dim>& solid_momentum) {
    nsolid_momentum_ += solid_momentum;
  }

  // assign nodal momentum from particles
  //! param[in] momentum momentum from particles : OK
  void assign_multimaterial_nodal_momenta(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& solid_momentum) {
    nsolid_momenta_.row(mat_id) += solid_momentum;
  }
  

  // assign nodal pore pressure from particles
  //! param[in] pressureXmass value of (pressure x mass) from particles : OK
  void assign_nodal_pore_pressure(const double &pressureXmass) {
    npressure_ += pressureXmass;
  }

  // assign nodal pore pressure from particles
  //! param[in] pressureXmass value of (pressure x mass) from particles : OK
  void assign_nodal_multimaterial_pore_pressures(const unsigned mat_id, const double &pressureXmass) {
    npressures_(mat_id) += pressureXmass;
  }

  // assign body forces from particles
  //! param[in] bodyforce body force for each phase from particles : OK
  void assign_body_force(const Eigen::Matrix<double,1,dim> &mixture_bodyforce,
			 const Eigen::Matrix<double,1,dim> &water_bodyforce){
    nmixture_body_force_ += mixture_bodyforce;
    nwater_body_force_ += water_bodyforce;
  }

  // assign multimaterial body forces from particles
  //! param[in] material id and body force for each phase from particles : OK
  void assign_multimaterial_body_forces(const unsigned mat_id, const Eigen::Matrix<double,1,dim> &mixture_bodyforce,
			 const Eigen::Matrix<double,1,dim> &water_bodyforce){
    nmixture_body_forces_.row(mat_id) += mixture_bodyforce;
    nwater_body_forces_.row(mat_id) += water_bodyforce;
  }

  // assign traction forces from particles
  //! param[in] tractionforce traction force for each phase from particles : OK
  void assign_traction_force(const Eigen::Matrix<double,1,dim> &mixture_tractionforce, const Eigen::Matrix<double,1,dim> &water_tractionforce){
    nmixture_trac_force_ += mixture_tractionforce;
    nwater_trac_force_ += water_tractionforce;
  }

  // assign traction forces from particles
  //! param[in] tractionforce traction force for each phase from particles : OK
  void assign_multimaterial_traction_forces(const unsigned mat_id, const Eigen::Matrix<double,1,dim> &mixture_tractionforce,
       const Eigen::Matrix<double,1,dim> &water_tractionforce){
    nmixture_trac_forces_.row(mat_id) += mixture_tractionforce;
    nwater_trac_forces_.row(mat_id) += water_tractionforce;
  }

  // assign solid internal force from particles
  //! param[in] intforce solid internal force from particles : OK
  void assign_internal_force(const Eigen::Matrix<double,1,dim> &mixture_intforce, const Eigen::Matrix<double,1,dim> &water_intforce){
    nmixture_int_force_ += mixture_intforce;
    nwater_int_force_ += water_intforce;
  }
  
  // assign solid internal force from particles
  //! param[in] intforce solid internal force from particles : OK
  void assign_multimaterial_internal_forces(const unsigned mat_id, const Eigen::Matrix<double,1,dim> &mixture_intforce,
       const Eigen::Matrix<double,1,dim> &water_intforce){
    nmixture_int_forces_.row(mat_id) += mixture_intforce;
    nwater_int_forces_.row(mat_id) += water_intforce;
  }
  
  void matrix_test();

  void write_para_data(std::ostream& oFile);
  // compute the damping forces
  void compute_nodal_damping_forces(const double &alpha) {};

  // compute nodal velocity from the momentum : OK
  void compute_nodal_initial_velocity(bool contact);
  // compute nodal pressure from the mapped values : OK
  void compute_nodal_pressure();

  // compute intermediate values
  void assign_nodal_intrmd_mass(const double& intrmd_mass) {
    nintrmd_mass_ += intrmd_mass;
  }
  void assign_multimaterial_nodal_intrmd_masses(const unsigned mat_id, const double& intrmd_mass) {
    nintrmd_masses_(mat_id) += intrmd_mass;
  }
  void compute_intermediate_solid_acceleration(const double& dt);
  void compute_intermediate_water_velocity();
  void compute_final_solid_acceleration(const double& dt);
  void compute_sp_final_solid_acceleration(const unsigned index, const double& dt);

  // compute multimaterial relative velocities
  void compute_multimaterial_relative_velocities();

  // set final pressure : End of time step
  void update_nodal_pressure_increment(const double &pressure) {
    npressure_increment_ = pressure;
    if (std::get<0>(npressure_constraints_)) 
      npressure_increment_ = 0.0;
  }

  // set final pressure : End of time step
  void update_multimaterial_nodal_pressure_increments(const unsigned mat_id, const double &pressure) {
    npressure_increments_(mat_id) = pressure;
    if (std::get<0>(npressure_constraints_)) 
      npressure_increments_(mat_id) = 0.0;
  }
 
  void assign_KS2_force(const Eigen::Matrix<double,1,dim>& force) {
    KS2_force_ = force;
  }

  void assign_multimaterial_KS2_forces(const unsigned mat_id, const Eigen::Matrix<double,1,dim>& force) {
    KS2_forces_.row(mat_id) = force;
  }

  //apply contact mechanics
  void apply_contact_mechanics(const double& dt);

  //moving mesh
  void update_mesh_configuration(double rigid_displacement, double soil_depth);

  // give id
  //! param[out] nid_ id of the node
  unsigned give_id() const {
      return nid_;
  }
  // give nodal nphase
  //! param[out] nphase_ number of phase of the node
  unsigned give_node_phase() const{
      return nphase_;
  }

  // give coordinates
  //! param[out] ncoord_ coordinates of the node
  Eigen::Matrix<double,1,dim> give_node_coordinates() const {
      return ncoord_;
  }

  // give penalty_factor
  //! param[out] penalty_factor_
  double give_node_penalty_factor() const {
    return penalty_factor_;
  }

  // give nodal solid mass
  //! param[out] nsolid_mass_ lumped solid mass of the node
  double give_node_solid_mass() const {
      return nsolid_mass_;
  }

  // give nodal water mass
  //! param[out] nwater_mass_ lumped water mass of the node
  double give_node_water_mass() const {
      return nwater_mass_;
  }

  // give nodal mixture mass
  //! param[out] nwater_mass_ lumped water mass of the node
  double give_node_mixture_mass() const {
      return nmixture_mass_;
  }

  // give node solid velocity
  //! param[out] nsolid_velocity_ solid velocity of the node
  Eigen::Matrix<double,1,dim> give_node_solid_initial_velocity() const {
      return nsolid_initial_velocity_;
  }

  Eigen::Matrix<double,1,dim> give_node_solid_multimaterial_initial_velocity(unsigned mat_id) const {
      return nsolid_initial_velocities_.row(mat_id);
  }

  Eigen::Matrix<double,1,dim> give_node_solid_final_velocity() const {
      return nsolid_final_velocity_;
  }

  Eigen::Matrix<double,1,dim> give_node_solid_increment_velocity() const {
    return (nsolid_final_velocity_ - nsolid_initial_velocity_);
  }
  
  Eigen::Matrix<double,1,dim> give_node_solid_multimaterial_increment_velocity(unsigned mat_id) const {
      return (nsolid_final_velocities_.row(mat_id) - nsolid_initial_velocities_.row(mat_id));
  }

  Eigen::Matrix<double,1,dim> give_node_solid_multimaterial_final_velocity(unsigned mat_id) const {
      return nsolid_final_velocities_.row(mat_id);
  }

  Eigen::Matrix<double,1,dim> give_node_solid_multimaterial_final_acceleration(unsigned mat_id) const {
      return nsolid_final_accelerations_.row(mat_id);
  }

  Eigen::Matrix<double,1,dim> give_node_solid_intermediate_velocity() const {
    return nsolid_int_velocity_;
  }
  Eigen::Matrix<double,1,dim> give_node_multimaterial_normal_vectors(unsigned mat_id) const{
    return normal_unit_vectors_.row(mat_id);
  }
  Eigen::Matrix<double,1,dim> give_node_multimaterial_solid_intermediate_velocities(unsigned mat_id) const {
    return nsolid_int_velocities_.row(mat_id);
  }

  Eigen::Matrix<double,1,dim> give_node_water_intermediate_velocity() const {
    return nwater_int_velocity_;
  }

  Eigen::Matrix<double,1,dim> give_node_multimaterial_water_intermediate_velocities(unsigned mat_id) const {
    return nwater_int_velocities_.row(mat_id);
  }

  // give node water velocity
  //! param[out] nwater_velocity_ water velocity of the node
  Eigen::Matrix<double,1,dim> give_node_water_initial_velocity() const {
      return nwater_initial_velocity_;
  }
  Eigen::Matrix<double,1,dim> give_node_water_final_velocity() const {
      return nwater_final_velocity_;
  }


  // give node solid acceleration
  //! param[out] nsolid_acceleration_ solid acceleration of the node
  Eigen::Matrix<double,1,dim> give_node_solid_final_acceleration() const {
    return (nsolid_final_acceleration_);
  }

  Eigen::Matrix<double,1,dim> give_node_multimaterial_solid_intermediate_accelerations(unsigned mat_id) const {
    return nsolid_int_accelerations_.row(mat_id);
  }

  // give node water acceleration
  //! param[out] nwater_acceleration_ water acceleration of the node
  Eigen::Matrix<double,1,dim> give_node_water_final_acceleration() const {
      return nwater_final_acceleration_;
  }

  // give nodal pore pressure
  //! param[out] npressure_ pore pressure of the node
  double give_node_pressure(unsigned mat_id) const {
      return npressures_(mat_id);
  }

  // give nodal pore pressure increment
  //! param[out] npressure_increment_ pore pressure increment of the node
  double give_node_pressure_increment() const {
      return npressure_increment_;
  }

  // give nodal pore pressure increment
  //! param[out] npressure_increment_ pore pressure increment of the node
  double give_node_multimaterial_pressure_increments(unsigned mat_id) const {
      return npressure_increments_(mat_id);
  }

  // give nodal multimaterial mixture mass
  //! param[out] mixture_masses_ multimaterial mixture mass
  double give_node_multimaterial_mixture_masses(unsigned mat_id) const {
      return nmixture_masses_(mat_id);
  }

  // give nodal solid body force
  //! param[out] nsolid_body_force_
  Eigen::Matrix<double,1,dim> give_node_mixture_body_force() const {
      return nmixture_body_force_;
  }

  // give nodal multimaterial mixture body force
  //! param[out] nmixture_body_forces
  Eigen::Matrix<double,1,dim> give_node_multimaterial_mixture_body_forces(unsigned mat_id) const {
    return nmixture_body_forces_.row(mat_id);
  }
  
  // give nodal multimaterial mixture body force
  //! param[out] nmixture_body_forces
  Eigen::Matrix<double,1,dim> give_node_multimaterial_mixture_trac_forces(unsigned mat_id) const {
    return nmixture_trac_forces_.row(mat_id);
  }

  // give nodal multimaterial internal body force
  //! param[out] nmixture_int_forces
  Eigen::Matrix<double,1,dim> give_node_multimaterial_mixture_internal_forces(unsigned mat_id) const {
    return nmixture_int_forces_.row(mat_id);
  }

  // give nodal water body force
  //! param[out] nwater_body_force_
  Eigen::Matrix<double,1,dim> give_node_water_body_force() const {
      return nwater_body_force_;
  }

  // give node solid traction force
  //! param[out] nsolid_trac__force_ solid traction force of the node
  Eigen::Matrix<double,1,dim> give_node_mixture_trac_force() const {
    return nmixture_trac_force_;
  }

  // give node water traction force
  //! param[out] nwater_trac__force_ solid traction force of the node
  Eigen::Matrix<double,1,dim> give_node_water_trac_force() const {
    return nwater_trac_force_;
  }

  // give node solid internal force
  //! param[out] nsolid_int_force_ solid internal force of the node
  Eigen::Matrix<double,1,dim> give_node_mixture_internal_force() const {
      return nmixture_int_force_;
  }

  // give node water internal force
  //! param[out] nwater_int_force_ water internal force of the node
  Eigen::Matrix<double,1,dim> give_node_water_internal_force() const {
      return nwater_int_force_;
  }

  // give nodal damping force fluid phase
  //! param[out] nwater_damping_
  Eigen::Matrix<double,1,dim> give_node_water_damping_force() {
    return nwater_damping_;
  }

  // give nodal damping force mixture
  //! param[out] nmixture_damping_
  Eigen::Matrix<double,1,dim> give_node_mixture_damping_force() {
    return nmixture_damping_;
  }

  // give node mixture forces
  //! param[out] mixture_force
  Eigen::Matrix<double,1,dim> give_node_mixture_forces() {    
    Eigen::Matrix<double,1,dim> mixture_force = nmixture_body_force_ + nmixture_int_force_ + nmixture_trac_force_;
    return mixture_force;
  }

  // give node fluid forces
  //! param[out] fluid_force
  Eigen::Matrix<double,1,dim> give_node_fluid_forces() {   
    Eigen::Matrix<double,1,dim> fluid_force = nwater_body_force_ + nwater_trac_force_ + nwater_int_force_;
    return fluid_force;
  }

  // give node KS2 force
  //! param[out] KS2_force
  Eigen::Matrix<double,1,dim> give_node_KS2_forces() {
    return KS2_force_;
  }

  // give solid velocity boundary condition 
  //! param[out] nsolid_constraints_
  std::vector<std::tuple<int,double> > give_node_solid_constraints() const {
      return nsolid_constraints_;
  }
  std::vector<std::tuple<int,double> > give_node_solid_accel_constraints() const {
      return nsolid_accel_constraints_;
  }

  // give water velocity boundary condition 
  //! param[out] nwater_constraints_
  std::vector<std::tuple<int,double> > give_node_water_constraints() const {
      return nwater_constraints_;
  }
  std::vector<std::tuple<int,double> > give_node_water_accel_constraints() const {
    return nwater_accel_constraints_;
  }

  // give pressure boundary condition
  //! param[out] the prescribed pressure
  std::tuple<bool,double> give_node_pressure_constraint() const {
      return npressure_constraints_;
  }

  // apply constraints to the solid velocity
  void apply_solid_velocity_constraint(Eigen::Matrix<double,1,dim> velocity);

  // apply constraints to the water velocity
  void apply_water_velocity_constraint(Eigen::Matrix<double,1,dim> velocity);

  // apply constraints to the solid acceleration
  void apply_solid_acceleration_constraint(Eigen::Matrix<double,1,dim> acceleration);

  // apply constraints to the water acceleration
  void apply_water_acceleration_constraint(Eigen::Matrix<double,1,dim> &acceleration);

  // apply constraints to the pressure
  void apply_pressure_constraint();

  // apply undrained conditions
  void apply_undrained_conditions(Eigen::Matrix<double,1,dim> water_relative_velocity);

  void set_moving_undrained_boundary() {
    std::get<0>(undrained_node_) = 1;
    std::get<1>(undrained_node_) = 1;
    moving_undrained_ = true;
    undrained_status_ = true;    
  }

  //! Material ids whose information was passed to this node
  std::set<unsigned> material_ids_;

private:
  // check the double precision of a value
  //! param[in] value 
  void check_double_precision(double& value);

public:
  // coordinates of node
  Eigen::Matrix<double,1,dim> ncoord_;

protected:
  // node ID
  unsigned nid_;

  // nodal solid velocity constraints
  std::vector<std::tuple<int, double> > nsolid_constraints_;
  std::vector<std::tuple<int, double> > nsolid_accel_constraints_;
  // nodal water velocity constraints
  std::vector<std::tuple<int, double> > nwater_constraints_;
  std::vector<std::tuple<int, double> > nwater_accel_constraints_;
  // nodal pressure constraints
  std::tuple<bool, double> npressure_constraints_;
  // undrained node: undrained direction; x=0 y =1 undrained normal +1 / -1
  std::tuple<unsigned,int> undrained_node_;
  bool undrained_status_;
  bool moving_undrained_;

  // Number of phase
  unsigned nphase_ = 1;

  // penalty factor for contact
  double penalty_factor_;

  // nodal lumped mass
  double nsolid_mass_;
  double nwater_mass_;
  double nmixture_mass_;
  double nintrmd_mass_;

  // multimaterial masses
  Eigen::Matrix<double,numMats,1> nsolid_masses_;
  Eigen::Matrix<double,numMats,1> nwater_masses_;
  Eigen::Matrix<double,numMats,1> nmixture_masses_;
  Eigen::Matrix<double,numMats,1> nintrmd_masses_;

  // nodal momentum
  Eigen::Matrix<double,1,dim> nsolid_momentum_;
  Eigen::Matrix<double,1,dim> nwater_momentum_;

  // multimaterial momenta
  Eigen::Matrix<double,numMats,dim> nsolid_momenta_;

  // nodal velocity
  Eigen::Matrix<double,1,dim> nsolid_initial_velocity_;
  Eigen::Matrix<double,1,dim> nwater_initial_velocity_;
  Eigen::Matrix<double,1,dim> nsolid_final_velocity_;
  Eigen::Matrix<double,1,dim> nwater_final_velocity_;

  // multimaterial velocities
  Eigen::Matrix<double,numMats,dim> nsolid_initial_velocities_;
  Eigen::Matrix<double,numMats,dim> nsolid_final_velocities_;
  
  // nodal acceleration
  Eigen::Matrix<double,1,dim> nsolid_final_acceleration_;
  Eigen::Matrix<double,1,dim> nwater_final_acceleration_;

  // multimaterial accelerations
  Eigen::Matrix<double,numMats,dim> nsolid_final_accelerations_;

  // nodal intermediate velocity
  Eigen::Matrix<double,1,dim> nsolid_int_acceleration_;
  Eigen::Matrix<double,1,dim> nsolid_int_velocity_;
  Eigen::Matrix<double,1,dim> nwater_int_velocity_;
  
  // multimaterial intermediate velocities
  Eigen::Matrix<double,numMats,dim> nsolid_int_accelerations_;
  Eigen::Matrix<double,numMats,dim> nsolid_int_velocities_;
  Eigen::Matrix<double,numMats,dim> nwater_int_velocities_;

  // nodal divergence free velocity
  Eigen::Matrix<double,1,dim> ndiv_free_velocity_;

  // multimaterial relative velocities
  Eigen::Matrix<double,numMats,dim> nsolid_relative_velocities_;

  // nodal pressure
  double npressure_;
  double npressure_increment_;

  // multimaterial pressures
  Eigen::Matrix<double,numMats,1> npressures_;
  Eigen::Matrix<double,numMats,1> npressure_increments_;

  // nodal body force
  Eigen::Matrix<double,1,dim> nmixture_body_force_;
  Eigen::Matrix<double,1,dim> nwater_body_force_;
  // multimaterial body forces
  Eigen::Matrix<double,numMats,dim> nmixture_body_forces_;
  Eigen::Matrix<double,numMats,dim> nwater_body_forces_;
  // nodal traction force
  Eigen::Matrix<double,1,dim> nmixture_trac_force_;
  Eigen::Matrix<double,1,dim> nwater_trac_force_;
  // multimaterial traction forces
  Eigen::Matrix<double,numMats,dim> nmixture_trac_forces_;
  Eigen::Matrix<double,numMats,dim> nwater_trac_forces_;
  // nodal internal force
  Eigen::Matrix<double,1,dim> nmixture_int_force_;
  Eigen::Matrix<double,1,dim> nwater_int_force_;
  // multimaterial internal forces
  Eigen::Matrix<double,numMats,dim> nmixture_int_forces_;
  Eigen::Matrix<double,numMats,dim> nwater_int_forces_;

  Eigen::Matrix<double,1,dim> nwater_damping_;
  Eigen::Matrix<double,1,dim> nmixture_damping_;

  Eigen::Matrix<double,1,dim> KS2_force_;
  Eigen::Matrix<double,numMats,dim> KS2_forces_;

  Eigen::Matrix<double,numMats,dim> domain_gradients_;
  Eigen::Matrix<double,numMats,dim> normal_unit_vectors_;

  bool assigned_; // this gives if the node was solved in projection solver
    // has to be improved

  //! TEMPORARY
  bool moving_boundary_;
  bool fixed_boundary_;

  bool free_pressure_;
  bool pressure_status_;

  unsigned step;
};

#include "Node.ipp"

#endif // MPM_NODE_H_
