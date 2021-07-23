mpm::material::MohrCoulomb::MohrCoulomb() {
    setProperty("density", density_);
    setProperty("porosity", porosity_);
    setProperty("permeability", permeability_);
    setProperty("youngModulus", E_);
    setProperty("poissonRatio", mu_);
    setProperty("frictionAngle", phi_);
    setProperty("cohesion", c_);
    setProperty("dilationAngle", psi_);
    setProperty("residualfrictionAngle", phi_resd_);
    setProperty("residualCohesion", c_resd_);
    setProperty("residualdilationAngle", psi_resd_);
    setProperty("tensileStrength", t_);
    setProperty("peakPlasticShearStrain", epds_peak_);
    setProperty("criticalPlasticShearStrain", epds_crit_);

    // initialise the plastic deviatoric strain vector
    PDS_ = Eigen::Matrix<double,6,1>::Zero();
    epds_ = 0.;
    phi_ = phi_ * (pi / 180.);
    phi_resd_ = phi_resd_ * (pi / 180.);
    psi_ = psi_ * (pi / 180.);
    psi_resd_ = psi_resd_ * (pi / 180.);
    this->compute_elastic_stiffness_matrix();
}

void mpm::material::MohrCoulomb::compute_elastic_stiffness_matrix() {
  double K, G;
  double alpha1, alpha2;
  De_ = Eigen::Matrix<double,6,6>::Zero();

  K = E_ / ( 3.0 * ( 1.0 - 2.0 * mu_ ) );
  G = E_ / ( 2.0 * (1.0 + mu_ ) );

  alpha1 = K + ( 4.0 * G / 3.0 );
  alpha2 = K - ( 2.0 * G / 3.0 );

  De_(0,0) = alpha1;   De_(0,1) = alpha2;   De_(0,2) = alpha2;
  De_(1,0) = alpha2;   De_(1,1) = alpha1;   De_(1,2) = alpha2;
  De_(2,0) = alpha2;   De_(2,1) = alpha2;   De_(2,2) = alpha1;
  De_(3,3) = G;
  De_(4,4) = G;
  De_(5,5) = G;

}

//! plane strain
//! Mohr Coulomb model is implemented based on the work by Dilan ROBERT
// param[in] dstrain: strain increment
// param[in] stress: stress
// param[in] epds: equivalent plastic deviatoric strain
// param[out] stress: stress is updated
// param[out] epds: equivalent plastic deviatoric strain is updated
void mpm::material::MohrCoulomb::compute_stress(const Eigen::Matrix<double,1,dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& epds) {

  VectorD6x1 dStrain = Eigen::Matrix<double,6,1>::Zero();  
  VectorD6x1 dStress = Eigen::Matrix<double,6,1>::Zero();
  VectorD6x1 totalStress = stress.transpose();
  VectorD6x1 devStress = Eigen::Matrix<double,6,1>::Zero();
  double j2, j3;
  double yieldF, yieldF_trial;
  bool yieldState;
  double pSoftening;
  double lambda, lambda_trial, pMultiplier;
  //epds_ = epds;

  // incremental strain
  dStrain(0) = dstrain(0); dStrain(1) = dstrain(1); dStrain(3) = dstrain(2);
  // softening: update the MC parameters
  if ((epds_peak_ - epds) >= 0.) {
    phi_cur_ = phi_; c_cur_ = c_; psi_cur_ = psi_;
  }
  else if ((epds - epds_peak_) > 0. || (epds_crit_- epds) >= 0.) {
    phi_cur_ = phi_resd_ + ((phi_ - phi_resd_) * (epds - epds_crit_) / (epds_peak_ - epds_crit_));
    psi_cur_ = psi_resd_ + ((psi_ - psi_resd_) * (epds - epds_crit_) / (epds_peak_ - epds_crit_));
    c_cur_ = c_resd_ + ((c_ - c_resd_) * (epds - epds_crit_) / (epds_peak_ - epds_crit_));
  }
  else if (epds > epds_crit_) {
    phi_cur_ = phi_resd_; psi_cur_ = psi_resd_; c_cur_ = c_resd_;
  }

  double meanP = (totalStress(0) + totalStress(1) + totalStress(2)) / 3.;
  if(meanP >= 0.0)
    meanP = 1.0;
  devStress(0) = totalStress(0) - meanP;
  devStress(1) = totalStress(1) - meanP;
  devStress(2) = totalStress(2) - meanP;
  devStress(3) = totalStress(3);
  devStress(4) = totalStress(4);
  devStress(5) = totalStress(5);
  if(dim == 2)
    devStress(4) = devStress(5) = 0.;

  // compute rho, theta and epsilon from the current stress state
  this->compute_rho_theta(totalStress, j2, j3);
  epsilon_ = (1. / sqrt(3)) * (stress(0) + stress(1) + stress(2));
  // compute yield function for the current stress state
  yieldF = sqrt(3. / 2.) * rho_ * ( (sin(theta_ + onethirdpi) / (sqrt(3.)*cos(phi_cur_))) + (cos(theta_ + onethirdpi) * tan(phi_cur_) / 3.) ) + (epsilon_ / 3.)*tan(phi_cur_) - c_cur_;
  if(yieldF > 1.E-15)
    yieldState = true;
  else
    yieldState = false;
  // compute plastic multiplier from the current stress state
  Eigen::Matrix<double,6,2> dFdP=this->compute_dFdP(j2,j3,devStress,pSoftening);
  Eigen::Matrix<double,1,6> dFdSigma = (dFdP.col(0)).transpose();
  Eigen::Matrix<double,1,6> dPdSigma = (dFdP.col(1)).transpose();
  lambda = ((dFdSigma * De_).dot(dStrain)) / ((dFdSigma*De_).dot(dPdSigma) + pSoftening);
  if(!yieldState)
    lambda = 0.;

  // compute the trial stress sigma_trial = sigma + De * dstrain
  VectorD6x1 trial_stress = stress.transpose() + (De_ * dStrain);
  meanP = (trial_stress(0) + trial_stress(1) + trial_stress(2)) / 3.;
  
  devStress(0) = trial_stress(0) - meanP;
  devStress(1) = trial_stress(1) - meanP;
  devStress(2) = trial_stress(2) - meanP;
  devStress(3) = trial_stress(3);
  devStress(4) = trial_stress(4);
  devStress(5) = trial_stress(5);
  if(dim == 2)
    devStress(4) = devStress(5) = 0.;

  this->compute_rho_theta(trial_stress, j2, j3);
  epsilon_ = (1. / sqrt(3)) * (trial_stress(0) + trial_stress(1) + trial_stress(2));
  yieldF_trial = sqrt(3. / 2.) * rho_ * ( (sin(theta_ + onethirdpi) / (sqrt(3.)*cos(phi_cur_))) + (cos(theta_ + onethirdpi) * tan(phi_cur_) / 3.) ) + (epsilon_ / 3.)*tan(phi_cur_) - c_cur_;

  bool yieldState_trial;
  if(yieldF_trial > 1.E-15)
    yieldState_trial = true;
  else
    yieldState_trial = false;
  Eigen::Matrix<double,6,2>  dFdP_trial=this->compute_dFdP(j2,j3,devStress,pSoftening);
  Eigen::Matrix<double,1,6> dFdSigma_trial = (dFdP_trial.col(0)).transpose();
  Eigen::Matrix<double,1,6>  dPdSigma_trial = (dFdP_trial.col(1)).transpose();

  lambda_trial = yieldF_trial / ((dFdSigma_trial*De_).dot(dPdSigma_trial) + pSoftening);

  if(yieldState)
    pMultiplier = lambda;
  else if (!yieldState) {
    if(yieldState_trial)
      pMultiplier = lambda_trial;
    else if (!yieldState_trial)
      pMultiplier = 0.;
  }

  // update stress (plastic correction)
  VectorD6x1 stress_new = trial_stress -  (pMultiplier * De_ * dPdSigma.transpose());
  stress = stress_new.transpose();
  // compute plastic deviatoric strain
  dStress = -totalStress + stress.transpose();
  VectorD6x1 dPstrain = dStrain - ((De_.inverse()) * dStress);
  if(dim ==2)
    dPstrain(4) = dPstrain(5) = 0.;
  PDS_ += dPstrain;

  // compute equivalent plastic deviatoric strain
  double epds_inc = (2./3.)*sqrt(0.5*((dPstrain(0)-dPstrain(1))*(dPstrain(0)-dPstrain(1)) + (dPstrain(1)-dPstrain(2))*(dPstrain(1)-dPstrain(2)) + (dPstrain(0)-dPstrain(2))*(dPstrain(0)-dPstrain(2))) + 3.*((dPstrain(3)*dPstrain(3)) + (dPstrain(4)*dPstrain(4)) + (dPstrain(5)*dPstrain(5))));
  epds_ += epds_inc;
  epds = epds_;
}


void mpm::material::MohrCoulomb::compute_rho_theta(const VectorD6x1 stress,
						   double& _j2, double& _j3) {
  double _thetaVal;
  double _meanP = (stress(0) + stress(1) + stress(2)) / 3.;
  VectorD6x1 _devstress = stress;
  _devstress(0) -= _meanP;
  _devstress(1) -= _meanP;
  _devstress(2) -= _meanP;
  if(dim ==2)
    _devstress(4) = _devstress(5) = 0.;
  
  // compute j2
  _j2 = (pow((stress(0) - stress(1)),2) + pow((stress(1) - stress(2)),2) + pow((stress(0) - stress(2)),2)) / 6.0 + pow(stress(3),2);
  if(dim == 3)
    _j2 += pow(stress(4),2) + pow(stress(5),2);
  // compute j3
  _j3 = (_devstress(0)*_devstress(1)*_devstress(2)) - (_devstress(2)*pow(_devstress(3),2));
  if(dim == 3)
    _j3 += ((2 * _devstress(3)*_devstress(4)*_devstress(5)) - (_devstress(0) * pow(_devstress(4),2) - _devstress(1)*pow(_devstress(5),2)));

  // compute theta value
  _thetaVal = 0.;
  if(fabs(_j2) > 1.E-12)
    _thetaVal = (3. * sqrt(3) / 2.) * (_j3 / pow(_j2,1.5));


  if(_thetaVal > 0.99) _thetaVal = 1.0;
  if(_thetaVal < -0.99) _thetaVal = -1.0;
  theta_ = (1./3.) * acos(_thetaVal);
  if(theta_ > onethirdpi) theta_ = onethirdpi;
  if(theta_ < 0.0) theta_ = 0.;
  
  rho_ = sqrt(2 * _j2);
}

Eigen::Matrix<double,6,2> mpm::material::MohrCoulomb::compute_dFdP(const double _j2, const double _j3, const VectorD6x1 _devstress, double& _psoftening){
  Eigen::Matrix<double,6,2> _dFanddP = Eigen::Matrix<double,6,2>::Zero();
  // compute dF / dEpsilon,  dF / dRho, dF / dTheta
  double _dFdEpsilon = tan(phi_cur_) / sqrt(3.);
  double _dFdRho = sqrt(3. / 2.) * ( (sin(theta_ + onethirdpi) / (sqrt(3.)*cos(phi_cur_))) + (cos(theta_ + onethirdpi) * tan(phi_cur_) / 3.) );
  double _dFdTheta = sqrt(3. / 2.) * rho_ * ( (cos(theta_ + onethirdpi) / (sqrt(3.)*cos(phi_cur_))) - (sin(theta_ + onethirdpi) * tan(phi_cur_) / 3.) );
  
  VectorD6x1 _dEpsilondSigma, _dRhodSigma, _dThetadSigma;
  _dEpsilondSigma = _dRhodSigma = _dThetadSigma = VectorD6x1::Zero();

  // compute dEpsilon / dSigma
  _dEpsilondSigma(0) = _dEpsilondSigma(1) = _dEpsilondSigma(2) = 1./sqrt(3.);
  
  // compute dRho / dSigma
  double multiplier = 1.;
  if (fabs(rho_) > 0.)
    multiplier = 1. / rho_;
  _dRhodSigma = multiplier * _devstress;
  if(dim == 2)
    _dRhodSigma(4) = _dRhodSigma(5) = 0.;
  
  // compute dTheta / dSigma
  double _r = 0.;
  if(fabs(_j2) > 1.E-12)
    _r = (3. * sqrt(3) / 2.) * (_j3 / pow(_j2,1.5));
  double devider = 1 - (_r*_r);
  if(devider <= 0.) devider = 0.001;
  double _dThetadR = -1 / (3. * sqrt(devider));
  double _dRdJ2 = (-9*sqrt(3.) / 4.) * _j3;
  if(fabs(_j2) > 1.E-12)
    _dRdJ2 = _dRdJ2 / pow(_j2, 2.5);
  double _dRdJ3 = 1.5 * sqrt(3.);
  if(fabs(_j2) > 1.E-12)
    _dRdJ3 = _dRdJ3 / pow(_j2, 1.5);

  VectorD6x1 _dJ2dSigma = _devstress;
  VectorD6x1 _dJ3dSigma = VectorD6x1::Zero();
  Eigen::Matrix<double,3,1> _dev1,_dev2,_dev3;
  _dev1(0) = _devstress(0); _dev1(1) = _devstress(3); _dev1(2) = _devstress(5);
  _dev2(0) = _devstress(3); _dev2(1) = _devstress(1); _dev2(2) = _devstress(4);
  _dev3(0) = _devstress(5); _dev3(1) = _devstress(4); _dev3(2) = _devstress(2);
  _dJ3dSigma(0) = _dev1.dot(_dev1) - (2. / 3.) * _j2;
  _dJ3dSigma(1) = _dev2.dot(_dev2) - (2. / 3.) * _j2;
  _dJ3dSigma(2) = _dev3.dot(_dev3) - (2. / 3.) * _j2;
  _dJ3dSigma(3) = _dev1.dot(_dev2);
  if(dim == 3) {
    _dJ3dSigma(4) = _dev2.dot(_dev3);
    _dJ3dSigma(5) = _dev1.dot(_dev3);
  }
  _dThetadSigma = _dThetadR * ((_dRdJ2 * _dJ2dSigma)+(_dRdJ3 * _dJ3dSigma));
  if(dim == 2)
    _dThetadSigma(4) = _dThetadSigma(5) = 0.;

  VectorD6x1 _dFdSigma = (_dFdEpsilon * _dEpsilondSigma) + (_dFdRho * _dRhodSigma) + (_dFdTheta * _dThetadSigma);
  if(dim == 2)
    _dFdSigma(4) = _dFdSigma(5) = 0.;

  // compute dP/dSigma
  double _R_mc = (3. - sin(phi_cur_)) / (6 * cos(phi_cur_));
  double _e = (3. - sin(phi_cur_)) / (3. + sin(phi_cur_));
  if((_e - 0.5) < 0.) _e = 0.501;
  if((_e - 1.) > 0.) _e = 1.0;
  double _sqpart = (4.*(1 - _e*_e)*pow(cos(theta_),2)) + (5*_e*_e) - (4.*_e);
  if(_sqpart < 0.) _sqpart = 0.00001;
  double _R_mw_den = (2*(1-_e*_e)*cos(theta_)) + ((2*_e - 1)*sqrt(_sqpart));
  if(fabs(_R_mw_den) < 1.E-18) _R_mw_den = 0.001;
  double _R_mw_num = (4*(1-_e*_e)*pow(cos(theta_),2)) + pow((2*_e - 1),2);
  double _R_mw = (_R_mw_num / _R_mw_den) * _R_mc;

  double _xi = 0.1;
  double _omega = pow((_xi*c_cur_*tan(psi_cur_)),2) + pow((_R_mw*sqrt(3./2.)*rho_),2);
  if(_omega < 1.E-22) _omega = 0.001;

  double _L = _R_mw_num;
  double _M = _R_mw_den;
  double _dLdTheta = -8 * (1-_e*_e)*cos(theta_)*sin(theta_);
  double _dMdTheta = (-2* (1-_e*_e)*sin(theta_)) + (0.5*(2*_e-1)*_dLdTheta)/sqrt(_sqpart);
  double _dRmwdTheta = ((_M*_dLdTheta)-(_L*_dMdTheta)) / (_M*_M);

  double _dPdEpsilon = tan(psi_cur_) / sqrt(3.);
  double _dPdRho = 3.*rho_*_R_mw*_R_mw / (2.*sqrt(_omega));
  double _dPdTheta = (3.*rho_*rho_ * _R_mw * _R_mc * _dRmwdTheta) / (2.*sqrt(_omega));

  VectorD6x1 _dPdSigma = (_dPdEpsilon*_dEpsilondSigma) + (_dPdRho*_dRhodSigma) + (_dPdTheta*_dThetadSigma);

  // compute softening part
  double _dPhidPstrain = 0.;
  double _dCdPstrain = 0.;
  if(epds_ > epds_peak_ && epds_ < epds_crit_) {
    _dPhidPstrain = (phi_resd_ - phi_) / (epds_crit_ - epds_peak_);
    _dCdPstrain = (c_resd_ - c_) / (epds_crit_ - epds_peak_);
  }
  double _dFdPhi = sqrt(3./2.)*rho_*((sin(phi_cur_)*sin(theta_ + onethirdpi) / (sqrt(3.)*cos(phi_cur_)*cos(phi_cur_))) + (cos(theta_ + onethirdpi)/(3.*cos(phi_cur_)*cos(phi_cur_))) ) + (epsilon_ / (sqrt(3.)*cos(phi_cur_)*cos(phi_cur_)));

  double _dFdC = -1.;
  _psoftening = (-1) * ((_dFdPhi*_dPhidPstrain) + (_dFdC*_dCdPstrain)) * _dPdRho;

  _dFanddP.col(0) = _dFdSigma;
  _dFanddP.col(1) = _dPdSigma;
  return _dFanddP;
}


