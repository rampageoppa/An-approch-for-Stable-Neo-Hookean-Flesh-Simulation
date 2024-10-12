#include <assemble_stiffness.h>
#include <iostream>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, Eigen::Ref<const Eigen::VectorXd> v0, 
                     double C, double D, bool stable) {


    K.setZero();
    K.resize(3 * q.rows(), 3 * q.rows());

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> triplets;

    for(int i = 0; i < T.rows(); i++)
    {

        Eigen::Matrix1212d H;
        d2V_linear_tetrahedron_dq2(H, q, V, T.row(i), v0(i), C, D, stable);


        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                for (int z = 0; z < 3; z++){
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k), H(3 * j + z, 3 * k)));
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k) + 1, H(3 * j + z, 3 * k + 1)));
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k) + 2, H(3 * j + z, 3 * k + 2)));
                }
            }
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());
    K = -1.0 * K;
}
