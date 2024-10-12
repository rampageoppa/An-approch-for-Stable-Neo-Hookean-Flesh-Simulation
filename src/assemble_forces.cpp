#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0,
                     double C, double D, bool stable) {


    f.resize(3 * q.rows());
    f.setZero();
    // https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
    // E: the spring connectivity matrix
    // need to consider gravity?

    for(int i = 0; i < T.rows(); i++)
    {
        Eigen::Vector12d dv;
        dV_linear_tetrahedron_dq(dv, q, V, T.row(i), v0(i), C, D, stable);

        for (int j = 0; j < 4; j++) {
            f(3 * T(i, j)) -= dv(3 * j);
            f(3 * T(i, j) + 1) -= dv(3 * j + 1);
            f(3 * T(i, j) + 2) -= dv(3 * j + 2);
        }
    }
};