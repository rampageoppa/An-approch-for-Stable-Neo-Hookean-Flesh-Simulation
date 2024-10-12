#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D, bool stable) {


    if(stable)
    {
        psi = ((2.0 * D) * pow((2.0 * C) / (2.0 * D) - F(0,0) * F(1,1) * F(2,2)+ F(0,0) * F(1,2) * F(2,1) + F(0,1) * F(1,0) * F(2,2) - F(0,1) * F(1,2) * F(2,0) - F(0,2) * F(1,0) * F(2,1) + F(0,2) * F(1,1) * F(2,0)+1.0,2.0))/2.0+((2.0 * C) *( F(0,0) * F(0,0) + F(0,1) * F(0,1) + F(0,2) * F(0,2) + F(1,0) * F(1,0) + F(1,1) * F(1,1) + F(1,2) * F(1,2) + F(2,0) * F(2,0) + F(2,1) * F(2,1) + F(2,2) * F(2,2) - 3.0)) / 2.0;
    }
    else
    {
        psi = C * ((F.transpose() * F).trace() / pow(F.determinant(), 2. / 3) - 3) + D * pow(F.determinant() - 1, 2);
    }

}