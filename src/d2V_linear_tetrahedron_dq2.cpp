#include <d2V_linear_tetrahedron_dq2.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <d2psi_neo_hookean_dq2.h>
#include <quadrature_single_point.h>
#include <iostream>

void d2V_linear_tetrahedron_dq2(Eigen::Matrix1212d &H, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D, bool stable) {

   auto neohookean_linear_tet = [&](Eigen::Matrix1212d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X){
        
       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);

       Eigen::MatrixXd B = Eigen::MatrixXd::Zero(9, 12);
       for(int i = 0; i < 3; i++){
           for(int j = 0; j < 4; j++){
               for(int z = 0; z < 3; z++){
                   B(i + 3 * z, 3 * j + z) = dphi(j, i);
               }
           }
       }

       Eigen::MatrixXd temp_x(3, 4);
       for(int i = 0; i < 4; i++){
           temp_x.block<3, 1>(0, i) = q.segment<3>(3 * element(i));
       }
       Eigen::Matrix3d F = temp_x * dphi;

       Eigen::Matrix99d temp_dp2;
       d2psi_neo_hookean_dF2(temp_dp2, F, C, D, stable);

       dV = B.transpose() * temp_dp2 * B;
    };

    quadrature_single_point(H, q, element, volume, neohookean_linear_tet);  
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix1212d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 12; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }

    H = Evec * DiagEval * Evec.transpose();

}
