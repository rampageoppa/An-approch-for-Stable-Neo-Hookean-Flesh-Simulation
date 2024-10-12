#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D, bool stable) {


    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {

        Eigen::MatrixXd temp_x(3, 4);
        for(int i = 0; i < 4; i++){
           temp_x.block<3, 1>(0, i) = q.segment<3>(3 * element(i));
        }

        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);
        psi_neo_hookean(e, temp_x * dphi, C, D, stable);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);  
    
}