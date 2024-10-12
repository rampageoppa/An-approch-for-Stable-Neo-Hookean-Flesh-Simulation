#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {

    Eigen::Vector3d temp_x0 = V.row(element(0));
    Eigen::Vector3d temp_x1 = V.row(element(1));
    Eigen::Vector3d temp_x2 = V.row(element(2));
    Eigen::Vector3d temp_x3 = V.row(element(3));
    
    Eigen::MatrixXd T(3, 3);
    T << temp_x1 - temp_x0, temp_x2 - temp_x0, temp_x3 - temp_x0;

    Eigen::Vector3d temp_phi;
    temp_phi = T.inverse() * (x - temp_x0);

    phi(0) = 1 - temp_phi(0) - temp_phi(1) - temp_phi(2);
    phi(1) = temp_phi(0);
    phi(2) = temp_phi(1);
    phi(3) = temp_phi(2);
}