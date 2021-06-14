mpm::material::ILE::ILE() {
    setProperty("density", density_);
    setProperty("porosity", porosity_);
    setProperty("permeability", permeability_);
    setProperty("youngModulus", E_);
    setProperty("poissonRatio", mu_);
    this->compute_elastic_stiffness_matrix();
}

void mpm::material::ILE::compute_elastic_stiffness_matrix() {
    double K, G;
    double a1, a2;
    K = E_ / (3.0 * (1. - 2. * mu_));
    G = E_ / (2.0 * (1. + mu_));
    a1 = K + (4.0 / 3.0) * G;
    a2 = K - (2.0 / 3.0) * G;
    De_(0,0)=a1;
    De_(0,1)=a2;
    De_(0,2)=a2;
    De_(0,3)=0;
    De_(0,4)=0;
    De_(0,5)=0;
    De_(1,0)=a2;
    De_(1,1)=a1;
    De_(1,2)=a2;
    De_(1,3)=0;
    De_(1,4)=0;
    De_(1,5)=0;
    De_(2,0)=a2;
    De_(2,1)=a2;
    De_(2,2)=a1;
    De_(2,3)=0;
    De_(2,4)=0;
    De_(2,5)=0;
    De_(3,0)= 0;
    De_(3,1)= 0;
    De_(3,2)= 0;
    De_(3,3)=G;
    De_(3,4)=0;
    De_(3,5)=0;
    De_(4,0)= 0;
    De_(4,1)= 0;
    De_(4,2)= 0;
    De_(4,3)=0;
    De_(4,4)=G;
    De_(4,5)=0;
    De_(5,0)= 0;
    De_(5,1)= 0;
    De_(5,2)= 0;
    De_(5,3)=0;
    De_(5,4)=0;
    De_(5,5)=G;
}


void mpm::material::ILE::compute_stress(const Eigen::Matrix<double, 1, dof>& dstrain, Eigen::Matrix<double,1,6>& stress, double& pressure) {

  Eigen::Matrix<double, 1, 6> incremental_strain = Eigen::Matrix<double, 1, 6>::Zero();
  Eigen::Matrix<double, 1, 6> incremental_stress = Eigen::Matrix<double, 1, 6>::Zero();
    if (dim == 2) {
        incremental_strain(0) = dstrain(0);
        incremental_strain(1) = dstrain(1);
        incremental_strain(3) = dstrain(2);
        incremental_stress = De_ * incremental_strain.transpose();
    }

    stress += incremental_stress;
}
