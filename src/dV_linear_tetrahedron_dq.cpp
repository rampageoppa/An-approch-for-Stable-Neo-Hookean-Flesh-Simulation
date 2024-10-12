#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D, bool stable) {


   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);

       Eigen::MatrixXd B(9, 12);
       B.setZero();

       for (int i = 0; i < 3; i++) {
           for (int j = 0; j < 4; j++) {
               for (int z = 0; z < 3; z++){
                   B(i + 3 * z, 3 * j + z) = dphi(j, i);
               }
           }
       }

        Eigen::MatrixXd temp_x(3, 4);
        for(int i = 0; i < 4; i++){
           temp_x.block<3, 1>(0, i) = q.segment<3>(3 * element(i));
        }
        Eigen::Matrix3d F = temp_x * dphi;

        Eigen::Vector9d temp_df;
        dpsi_neo_hookean_dF(temp_df, F, C, D, stable);

        dV = B.transpose() * temp_df;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
    
}