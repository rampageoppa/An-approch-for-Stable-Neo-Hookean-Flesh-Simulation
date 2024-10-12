#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  F - the dense 3x3 deformation gradient
//  C,D - material parameters for the Neo-Hookean model
// stable - whether we are presenting the stable noe-hookean or not
//Output:
//  psi - the neohookean energy
void psi_neo_hookean(double &psi, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D, bool stable);