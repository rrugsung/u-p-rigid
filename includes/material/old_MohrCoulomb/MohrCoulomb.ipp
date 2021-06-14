mpm::material::MohrCoulomb::MohrCoulomb() {
  setProperty("density", density_);
  setProperty("porosity", porosity_);
  setProperty("permeability", permeability_);
  setProperty("youngModulus", E_);
  setProperty("poissonRatio", mu_);
  setProperty("frictionAngle", phi_);
  setProperty("cohesion", c_);
  setProperty("dilationAngle", psi_);
  setProperty("tensileStrength", t_);

  pstrain_.setZero();
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
//! Mohr Coulomb model is same as FLAC
void mpm::material::MohrCoulomb::compute_stress(const Eigen::Matrix<double, 1, dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& plastic_strain) {

  Eigen::Matrix<double,6,1> dStrain = Eigen::Matrix<double,6,1>::Zero();  
  Eigen::Matrix<double,6,1> dStress = Eigen::Matrix<double,6,1>::Zero();
  Eigen::Matrix<double,6,1> total_stress = Eigen::Matrix<double,6,1>::Zero();
  Eigen::Matrix<double,6,1> unity; unity << 1,1,1,0,0,0;

  double phi, psi;
  double Nphi, Npsi;
  double sigma1, sigma2, sigma3;
  double f_s, f_t;
  double lambda_s, lambda_t;
    lambda_s = 0;
  double dplastic_strain_1, dplastic_strain_2, dplastic_strain_3;

  phi_crit_ = 20.;
  c_crit_ = 50;
  gamma_peak_ = 0.15;
  gamma_crit_ = 0.25;

  phi = ( phi_ * pi ) / 180;
  psi = ( psi_ * pi ) / 180;

  Nphi = ( 1.0 + sin(phi) ) / ( 1.0 - sin(phi) );
  Npsi = ( 1.0 + sin(psi) ) / ( 1.0 - sin(psi) );

  //! initial guess of elastic stress 
  this->compute_elastic_stiffness_matrix();
  dStrain(0) = dstrain(0);       dStrain(1) = dstrain(1);      dStrain(3) = dstrain(2);
  dStress = De_ * dStrain;
  total_stress = stress.transpose() + dStress;// - (pressure * unity);

  pstrain_ += dstrain;
  double vol = (pstrain_(0) + pstrain_(1))/3.;
  double dev_1 = pstrain_(0) - vol;
  double dev_2 = pstrain_(1) - vol;
  double dev_3 = pstrain_(2);
  
  //plastic_strain = sqrt(0.5*( pow(dev_1,2) + pow(dev_2,2) + pow(dev_3,2) ) );

  if (plastic_strain > gamma_peak_ && plastic_strain < gamma_crit_) {
    if(phi_ < 40.)
      phi_ = phi_crit_ + ((phi_ - phi_crit_) * (gamma_crit_ - plastic_strain) / (gamma_crit_ - gamma_peak_));
   // if (c_ > 0.)
   //    c_ = c_crit_ + ((c_ - c_crit_) * (gamma_crit_ - plastic_strain) / (gamma_crit_ - gamma_peak_));
  }
  if (plastic_strain > gamma_crit_) {
    if(phi_ < 40.)
      phi_ = phi_crit_;
    // if(c_ > 0.0)
    // c_ = c_crit_;
  }

  double sigma_xx = total_stress(0);
  double sigma_yy = total_stress(1);
  double sigma_zz = total_stress(2);
  double sigma_xy = total_stress(3);
  
  //! mohr's circle centre / radius / angle
  double mohr_radius = std::sqrt( ( sigma_xy * sigma_xy ) + ( (sigma_xx - sigma_yy) * (sigma_xx - sigma_yy) / 4.0 ) );
  double mohr_centre = ( sigma_xx + sigma_yy ) / 2.0;

  //! maximum and minimum stress for 2D case
  double max_sigma = mohr_centre - mohr_radius;
  double min_sigma = mohr_centre + mohr_radius;
  unsigned stress_order;

  double sin_2theta, cos_2theta;
  cos_2theta = ( sigma_xx - sigma_yy ) / ( min_sigma - max_sigma );
  sin_2theta = ( 2.0 * sigma_xy ) / (min_sigma - max_sigma );

  if (sigma_zz > min_sigma) {
    sigma1 = max_sigma;
    sigma2 = min_sigma;
    sigma3 = sigma_zz;
    stress_order = 1;
  }
  else if (sigma_zz < max_sigma) {
    sigma1 = sigma_zz;
    sigma2 = max_sigma;
    sigma3 = min_sigma;
    stress_order = 2;
  }
  else {
    sigma1 = max_sigma;
    sigma2 = sigma_zz;
    sigma3 = min_sigma;
    stress_order = 3;
  }

  //! Mohr Coulomb yield function
  f_s = sigma1 - ( sigma3 * Nphi ) + ( 2.0 * c_ * std::sqrt(Nphi) );
  //! Tension yield function
  f_t = t_ - sigma3;

  //! h-function which defines the mode of failure (shear / tension)
  double alpha_p = std::sqrt(1.0 + Nphi * Nphi) + Nphi;
  double sigma_p = ( t_ * Nphi ) - ( 2.0 * c_ * std::sqrt(Nphi) );
  double h_function = sigma3 - t_ + ( alpha_p * ( sigma1 - sigma_p ) );

  double K = E_ / ( 3.0 * ( 1.0 - 2.0 * mu_ ) );
  double G = E_ / ( 2.0 * (1.0 + mu_ ) );
  double alpha1 = K + ( 4.0 * G / 3.0 );
  double alpha2 = K - ( 2.0 * G / 3.0 );

  if ( h_function < 0. ) { // shear failure
    if ( f_s < 0. ) {
      lambda_s = f_s / ( (alpha1 - alpha2 * Npsi ) - ( alpha2 - alpha1 * Npsi ) * Nphi );
      sigma1 -= ( lambda_s * ( alpha1 - alpha2 * Npsi ) );
      sigma2 -= ( lambda_s * alpha2 * ( 1.0 - Npsi ) );
      sigma3 -= ( lambda_s * ( alpha2 - alpha1 * Npsi ) );
      //std::cout << "shear failure" << "\n";
    }
  } else if (h_function > 0. ) { // tension failure
    if (f_t < 0. ) {
      lambda_t = f_t / alpha1;
      sigma1 += ( lambda_t * alpha2 );
      sigma2 += ( lambda_t * alpha2 );
      sigma3 += ( lambda_t * alpha1 );
    }
  }

  if (sigma1 == sigma3) {
    cos_2theta = 1.0;
    sin_2theta = 0.0;  
  }

  stress = Eigen::Matrix<double,1,6>::Zero();
  if (stress_order == 3) {
    mohr_centre = 0.5 * ( sigma1 + sigma3 );
    mohr_radius = 0.5 * ( sigma3 - sigma1 );
    stress(0) = mohr_centre + ( mohr_radius * cos_2theta );
    stress(1) = mohr_centre - ( mohr_radius * cos_2theta );
    stress(2) = sigma2;
    stress(3) = mohr_radius * sin_2theta;
  } else if (stress_order == 2) {
    mohr_centre = 0.5 * ( sigma2 + sigma3 );
    mohr_radius = 0.5 * ( sigma3 - sigma2 );
    stress(0) = mohr_centre + ( mohr_radius * cos_2theta );
    stress(1) = mohr_centre - ( mohr_radius * cos_2theta );
    stress(2) = sigma1;
    stress(3) = mohr_radius * sin_2theta; 
  } else if (stress_order == 1) {
    mohr_centre = 0.5 * ( sigma1 + sigma2 );
    mohr_radius = 0.5 * ( sigma2 - sigma1 );
    stress(0) = mohr_centre + ( mohr_radius * cos_2theta );
    stress(1) = mohr_centre - ( mohr_radius * cos_2theta );
    stress(2) = sigma3;
    stress(3) = mohr_radius * sin_2theta; 
  }
  dplastic_strain_1 = lambda_s;
  dplastic_strain_2 = 0;
  dplastic_strain_3 = -lambda_s * Npsi;
  double principal_inv1 = (sigma1 + sigma2+sigma3)/3.;
  double cauchy_inv1 = (stress(0)+stress(1)+stress(2))/3.;
  if(std::fabs(principal_inv1 - cauchy_inv1)>1.E-8)
    std::cout << "ERROR" << "\n";
  

  //plastic_strain += (sqrt(2) /3.) * sqrt(pow(dplastic_strain_1-dplastic_strain_2,2) + pow(dplastic_strain_2-dplastic_strain_3,2) + pow(dplastic_strain_3 - dplastic_strain_1,2));
  //stress += (pressure * unity);
  //stress += (pressure * unity);
}
