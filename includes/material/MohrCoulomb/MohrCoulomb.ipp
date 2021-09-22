mpm::material::MohrCoulomb::MohrCoulomb() {
    setProperty("density", density_);
    setProperty("porosity", porosity_);
    setProperty("permeability", permeability_);
    setProperty("youngModulus", youngs_modulus_);
    setProperty("poissonRatio", poisson_ratio_);
    setProperty("frictionAngle", phi_peak_);
    setProperty("cohesion", cohesion_peak_);
    setProperty("dilationAngle", psi_peak_);
    setProperty("residualfrictionAngle", phi_residual_);
    setProperty("residualCohesion", cohesion_residual_);
    setProperty("residualdilationAngle", psi_residual_);
    setProperty("tensileStrength", tension_cutoff_);
    setProperty("peakPlasticShearStrain", pdstrain_peak_);
    setProperty("criticalPlasticShearStrain", pdstrain_residual_);

    // initialise the plastic deviatoric strain vector
    PDS_ = Eigen::Matrix<double,6,1>::Zero();
    epds_ = 0.;
    phi_peak_ = phi_peak_ * (pi / 180.);
    phi_residual_ = phi_residual_ * (pi / 180.);
    psi_peak_ = psi_peak_ * (pi / 180.);
    psi_residual_ = psi_residual_ * (pi / 180.);

    // Bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    // Shear modulus
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));
    this->compute_elastic_stiffness_matrix();
}

void mpm::material::MohrCoulomb::compute_elastic_stiffness_matrix() {
  // Shear modulus
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  // compute elastic stiffness matrix
  // clang-format off
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;

}

//! plane strain
//! Mohr Coulomb model is implemented based on the work by Dilan ROBERT
// param[in] dstrain: strain increment
// param[in] stress: stress
// param[in] epds: equivalent plastic deviatoric strain
// param[out] stress: stress is updated
// param[out] epds: equivalent plastic deviatoric strain is updated
void mpm::material::MohrCoulomb::compute_stress(const Eigen::Matrix<double,1,dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& epds) {
  // Get equivalent plastic deviatoric strain
  const double pdstrain = epds;
  VectorD6x1 dStrain = Eigen::Matrix<double,6,1>::Zero();  
  VectorD6x1 TotalStress = stress.transpose();
  dStrain(0) = dstrain(0);
  dStrain(1) = dstrain(1);
  dStrain(3) = dstrain(2);

  // Update MC parameters using a linear softening rule
  if (pdstrain > pdstrain_peak_) {
    if (pdstrain < pdstrain_residual_) {
      phi_cur =
          phi_residual_ +
          ((phi_peak_ - phi_residual_) * (pdstrain - pdstrain_residual_) /
           (pdstrain_peak_ - pdstrain_residual_));
      psi_cur =
          psi_residual_ +
          ((psi_peak_ - psi_residual_) * (pdstrain - pdstrain_residual_) /
           (pdstrain_peak_ - pdstrain_residual_));
      cohesion_cur =
          cohesion_residual_ + ((cohesion_peak_ - cohesion_residual_) *
                                (pdstrain - pdstrain_residual_) /
                                (pdstrain_peak_ - pdstrain_residual_));
    } else {
      phi_cur = phi_residual_;
      psi_cur = psi_residual_;
      cohesion_cur = cohesion_residual_;
    }
  }
  else {
      phi_cur = phi_peak_;
      psi_cur = psi_peak_;
      cohesion_cur = cohesion_peak_;
  }
  //-------------------------------------------------------------------------
  // Elastic-predictor stage: compute the trial stress
  VectorD6x1 trial_stress = stress.transpose() + (this->de_ * dStrain);
  // Compute stress invariants based on trial stress
  this->compute_stress_invariants(trial_stress);
  // Compute yield function based on the trial stress
  Eigen::Matrix<double, 2, 1> yield_function_trial;
  auto yield_type_trial =
      this->compute_yield_state(&yield_function_trial, epsilon_, rho_, theta_, phi_cur, cohesion_cur);

  // Return the updated stress in elastic state
  if (yield_type_trial == mpm::material::MohrCoulomb::FailureState::Elastic)
    stress = trial_stress.transpose();
  
  else {
    //-------------------------------------------------------------------------
    // Plastic-corrector stage: correct the stress back to the yield surface
    // Define tolerance of yield function
    const double Tolerance = 1E-1;
    // Compute plastic multiplier based on trial stress (Lambda trial)
    double softening_trial = 0.;
    double dp_dq_trial = 0.;
    VectorD6x1 df_dsigma_trial = VectorD6x1::Zero();
    VectorD6x1 dp_dsigma_trial = VectorD6x1::Zero();
    this->compute_df_dp(yield_type_trial, trial_stress,
                        &df_dsigma_trial, &dp_dsigma_trial, rho_, theta_, phi_cur, psi_cur, pdstrain, &dp_dq_trial,
                        &softening_trial);
    double yield_trial = 0.;
    if (yield_type_trial == mpm::material::MohrCoulomb::FailureState::Tensile)
      yield_trial = yield_function_trial(0);
    if (yield_type_trial == mpm::material::MohrCoulomb::FailureState::Shear)
      yield_trial = yield_function_trial(1);
    double lambda_trial =
        yield_trial /
        ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
        softening_trial);
    // Compute stress invariants based on stress input
    this->compute_stress_invariants(TotalStress);
    // Compute yield function based on stress input
    Eigen::Matrix<double, 2, 1> yield_function;
    auto yield_type = this->compute_yield_state(&yield_function, epsilon_, rho_, theta_, phi_cur, cohesion_cur);

    // Initialise value of yield function based on stress
    double yield{std::numeric_limits<double>::max()};
    if (yield_type == mpm::material::MohrCoulomb::FailureState::Tensile)
      yield = yield_function(0);
    if (yield_type == mpm::material::MohrCoulomb::FailureState::Shear)
      yield = yield_function(1);
    // Compute plastic multiplier based on stress input (Lambda)
    double softening = 0.;
    double dp_dq = 0.;
    VectorD6x1 df_dsigma = VectorD6x1::Zero();
    VectorD6x1 dp_dsigma = VectorD6x1::Zero();
    this->compute_df_dp(yield_type, TotalStress, &df_dsigma, &dp_dsigma, rho_, theta_, phi_cur, psi_cur, pdstrain,
                        &dp_dq, &softening);
    const double lambda =
        ((df_dsigma.transpose() * this->de_).dot(dStrain)) /
        (((df_dsigma.transpose() * this->de_).dot(dp_dsigma)) + softening);
    // Initialise updated stress
    VectorD6x1 updated_stress = trial_stress;
    // Initialise incremental of plastic deviatoric strain
    double dpdstrain = 0.;
    // Correction stress based on stress
    if (fabs(yield) < Tolerance) {
      // Compute updated stress
      updated_stress -= (lambda * this->de_ * dp_dsigma);
      // Compute incremental of plastic deviatoric strain
      dpdstrain = lambda * dp_dq;
    } else {
      // Compute updated stress
      updated_stress -= (lambda_trial * this->de_ * dp_dsigma_trial);
      // Compute incremental of plastic deviatoric strain
      dpdstrain = lambda_trial * dp_dq_trial;
    }

    // Define the maximum iteration step
    const int itr_max = 100;
    // Correct the stress again
    for (unsigned itr = 0; itr < itr_max; ++itr) {
      // Check the update stress
      // Compute stress invariants based on updated stress
      this->compute_stress_invariants(updated_stress);
      // Compute yield function based on updated stress
      yield_type_trial =
          this->compute_yield_state(&yield_function_trial, epsilon_, rho_, theta_, phi_cur, cohesion_cur);
      // Check yield function
      if (yield_function_trial(0) < Tolerance &&
          yield_function_trial(1) < Tolerance) {
        break;
      }
      // Compute plastic multiplier based on updated stress
      this->compute_df_dp(yield_type_trial, updated_stress, &df_dsigma_trial, 
                          &dp_dsigma_trial, rho_, theta_, phi_cur, psi_cur, 
                          pdstrain, &dp_dq_trial, &softening_trial);
      if (yield_type_trial == mpm::material::MohrCoulomb::FailureState::Tensile)
        yield_trial = yield_function_trial(0);
      if (yield_type_trial == mpm::material::MohrCoulomb::FailureState::Shear)
        yield_trial = yield_function_trial(1);
      // Compute plastic multiplier based on updated stress
      lambda_trial =
          yield_trial /
          ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
           softening_trial);
      // Correct stress back to the yield surface
      updated_stress -= (lambda_trial * this->de_ * dp_dsigma_trial);
      // Update incremental of plastic deviatoric strain
      dpdstrain += lambda_trial * dp_dq_trial;
    }
    // Compute stress invariants based on updated stress
    this->compute_stress_invariants(updated_stress);
    // Update plastic deviatoric strain
    epds += dpdstrain;
    epds_ += epds;

    stress =  updated_stress.transpose();
  }
}

typename mpm::material::MohrCoulomb::FailureState
  mpm::material::MohrCoulomb::compute_yield_state(Eigen::Matrix<double, 2, 1>* yield_function, double epsilon, double rho, double theta, double phi, double cohesion){
  // Tolerance for yield function
  const double Tolerance = -1E-1;
  // Compute yield functions (tension & shear)
  // Tension
  (*yield_function)(0) = std::sqrt(2. / 3.) * cos(theta) * rho +
                         epsilon / std::sqrt(3.) - tension_cutoff_;
  // Shear
  (*yield_function)(1) =
      std::sqrt(1.5) * rho *
          ((sin(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) +
           (cos(theta + M_PI / 3.) * tan(phi) / 3.)) +
      (epsilon / std::sqrt(3.)) * tan(phi) - cohesion;
  // Initialise yield status (0: elastic, 1: tension failure, 2: shear failure)
  auto yield_type = mpm::material::MohrCoulomb::FailureState::Elastic;
  // Check for tension and shear
  if ((*yield_function)(0) > Tolerance && (*yield_function)(1) > Tolerance) {
    // Compute tension and shear edge parameters
    const double n_phi = (1. + sin(phi)) / (1. - sin(phi));
    const double sigma_p =
        tension_cutoff_ * n_phi - 2. * cohesion * std::sqrt(n_phi);
    const double alpha_p = std::sqrt(1. + n_phi * n_phi) + n_phi;
    // Compute the shear-tension edge
    const double h =
        (*yield_function)(0) +
        alpha_p * (std::sqrt(2. / 3.) * cos(theta - 4. * M_PI / 3.) * rho +
                   epsilon / std::sqrt(3.) - sigma_p);
    // Tension
    if (h > std::numeric_limits<double>::epsilon())
      yield_type = mpm::material::MohrCoulomb::FailureState::Tensile;
    // Shear
    else
      yield_type = mpm::material::MohrCoulomb::FailureState::Shear;
  }
  // Shear failure
  if ((*yield_function)(0) < Tolerance && (*yield_function)(1) > Tolerance)
    yield_type = mpm::material::MohrCoulomb::FailureState::Shear;
  // Tension failure
  if ((*yield_function)(0) > Tolerance && (*yield_function)(1) < Tolerance)
    yield_type = mpm::material::MohrCoulomb::FailureState::Tensile;

  return yield_type;
}
void mpm::material::MohrCoulomb::compute_stress_invariants(const VectorD6x1 stress) {
  // Compute the mean pressure
  epsilon_ = mpm::material::p(stress) * std::sqrt(3.);
  // Compute theta value
  theta_ = mpm::material::lode_angle(stress);
  // Compute rho
  rho_ = std::sqrt(2. * mpm::material::j2(stress));
}

void mpm::material::MohrCoulomb::compute_df_dp(mpm::material::MohrCoulomb::FailureState yield_type,
    const VectorD6x1& stress, VectorD6x1* df_dsigma, VectorD6x1* dp_dsigma, const double rho, const double theta, 
    const double phi, const double psi, double pdstrain, double* dp_dq, double* softening){

  // Compute dF / dEpsilon,  dF / dRho, dF / dTheta
  double df_depsilon, df_drho, df_dtheta;
  // Values in tension yield
  if (yield_type == mpm::material::MohrCoulomb::FailureState::Tensile) {
    df_depsilon = 1. / std::sqrt(3.);
    df_drho = std::sqrt(2. / 3.) * cos(theta);
    df_dtheta = -std::sqrt(2. / 3.) * rho * sin(theta);
  }
  // Values in shear yield / elastic
  else {
    df_depsilon = tan(phi) / std::sqrt(3.);
    df_drho = std::sqrt(1.5) *
              ((sin(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) +
               (cos(theta + M_PI / 3.) * tan(phi) / 3.));
    df_dtheta = std::sqrt(1.5) * rho *
                ((cos(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) -
                 (sin(theta + M_PI / 3.) * tan(phi) / 3.));
  }
  // Compute dEpsilon / dSigma
  VectorD6x1 depsilon_dsigma = mpm::material::dp_dsigma(stress) * std::sqrt(3.);
  // Initialise dRho / dSigma
  VectorD6x1 drho_dsigma = mpm::material::dq_dsigma(stress) * std::sqrt(2. / 3.);
  // Compute dtheta / dsigma
  VectorD6x1 dtheta_dsigma = mpm::material::dtheta_dsigma(
      stress, std::numeric_limits<double>::epsilon());
  // Compute dF/dSigma
  (*df_dsigma) = (df_depsilon * depsilon_dsigma) + (df_drho * drho_dsigma) +
                 (df_dtheta * dtheta_dsigma);
  // Compute dp/dsigma and dp/dj in tension yield
  if (yield_type == mpm::material::MohrCoulomb::FailureState::Tensile) {
    // Define deviatoric eccentricity
    const double et_value = 0.6;
    // Define meridional eccentricity
    const double xit = 0.1;
    // Compute Rt
    double sqpart = 4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                    5. * et_value * et_value - 4. * et_value;
    if (sqpart < std::numeric_limits<double>::epsilon()) sqpart = 1.E-5;
    double rt_den = 2. * (1 - et_value * et_value) * cos(theta) +
                    (2. * et_value - 1) * std::sqrt(sqpart);
    const double rt_num =
        4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
        (2. * et_value - 1) * (2. * et_value - 1);
    if (fabs(rt_den) < std::numeric_limits<double>::epsilon()) rt_den = 1.E-5;
    const double rt = rt_num / (3. * rt_den);
    // Compute dP/dRt
    const double dp_drt =
        1.5 * rho * rho * rt /
        std::sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                  1.5 * rt * rt * rho * rho);
    // Compute dP/dRho
    const double dp_drho =
        1.5 * rho * rt * rt /
        std::sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                  1.5 * rt * rt * rho * rho);
    // Compute dP/dEpsilon
    const double dp_depsilon = 1. / std::sqrt(3.);
    // Compute dRt/dThera
    const double drtden_dtheta =
        -2. * (1 - et_value * et_value) * sin(theta) -
        (2. * et_value - 1) * 4. * (1 - et_value * et_value) * cos(theta) *
            sin(theta) /
            std::sqrt(4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                      5. * et_value * et_value - 4. * et_value);
    const double drtnum_dtheta =
        -8. * (1 - et_value * et_value) * cos(theta) * sin(theta);
    const double drt_dtheta =
        (drtnum_dtheta * rt_den - drtden_dtheta * rt_num) /
        (3. * rt_den * rt_den);
    // Compute dP/dSigma
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_drt * drt_dtheta * dtheta_dsigma);
    // Compute dP/dJ
    (*dp_dq) = dp_drho * std::sqrt(2. / 3.);
  }
  // Compute dp/dsigma and dp/dj in shear yield
  else {
    // Compute Rmc
    const double r_mc = (3. - sin(phi)) / (6 * cos(phi));
    // Compute deviatoric eccentricity
    double e_val = (3. - sin(phi)) / (3. + sin(phi));
    if (e_val <= 0.5) e_val = 0.5 + 1.E-10;
    if (e_val > 1.) e_val = 1.;
    // Compute Rmw
    double sqpart = (4. * (1 - e_val * e_val) * std::pow(cos(theta), 2)) +
                    (5 * e_val * e_val) - (4. * e_val);
    if (sqpart < std::numeric_limits<double>::epsilon()) sqpart = 1.E-5;
    double m = (2. * (1 - e_val * e_val) * cos(theta)) +
               ((2. * e_val - 1) * std::sqrt(sqpart));
    if (fabs(m) < std::numeric_limits<double>::epsilon()) m = 1.E-5;
    const double l = (4. * (1. - e_val * e_val) * std::pow(cos(theta), 2)) +
                     std::pow((2. * e_val - 1.), 2);
    const double r_mw = (l / m) * r_mc;
    // Initialise meridional eccentricity
    const double xi = 0.1;
    double omega = std::pow((xi * cohesion_peak_ * tan(psi)), 2) +
                   std::pow((r_mw * std::sqrt(1.5) * rho), 2);
    if (omega < std::numeric_limits<double>::epsilon()) omega = 1.E-5;
    const double dl_dtheta =
        -8. * (1. - e_val * e_val) * cos(theta) * sin(theta);
    const double dm_dtheta =
        (-2. * (1. - e_val * e_val) * sin(theta)) +
        (0.5 * (2. * e_val - 1.) * dl_dtheta) / std::sqrt(sqpart);
    const double drmw_dtheta = ((m * dl_dtheta) - (l * dm_dtheta)) / (m * m);
    const double dp_depsilon = tan(psi) / std::sqrt(3.);
    const double dp_drho = 3. * rho * r_mw * r_mw / (2. * std::sqrt(omega));
    const double dp_dtheta =
        (3. * rho * rho * r_mw * r_mc * drmw_dtheta) / (2. * std::sqrt(omega));
    // compute the value of dp/dsigma and dp/dj in shear yield
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_dtheta * dtheta_dsigma);
    (*dp_dq) = dp_drho * std::sqrt(2. / 3.);
  }
  // Compute softening part
  double dphi_dpstrain = 0.;
  double dc_dpstrain = 0.;
  (*softening) = 0.;
  if (pdstrain > pdstrain_peak_ &&
      pdstrain < pdstrain_residual_) {
    // Compute dPhi/dPstrain
    dphi_dpstrain =
        (phi_residual_ - phi_peak_) / (pdstrain_residual_ - pdstrain_peak_);
    // Compute dc/dPstrain
    dc_dpstrain = (cohesion_residual_ - cohesion_peak_) /
                  (pdstrain_residual_ - pdstrain_peak_);
    // Compute dF/dPstrain
    double df_dphi =
        std::sqrt(1.5) * rho *
            ((sin(phi) * sin(theta + M_PI / 3.) /
              (std::sqrt(3.) * cos(phi) * cos(phi))) +
             (cos(theta + M_PI / 3.) / (3. * cos(phi) * cos(phi)))) +
        (mpm::material::p(stress) / (cos(phi) * cos(phi)));
    double df_dc = -1.;
    (*softening) =
        (-1.) * ((df_dphi * dphi_dpstrain) + (df_dc * dc_dpstrain)) * (*dp_dq);
  }
}


