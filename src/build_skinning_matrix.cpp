#include <build_skinning_matrix.h>
#include <phi_linear_tetrahedron.h>
#include <vector>
#include <iostream>

void build_skinning_matrix(Eigen::SparseMatrixd &N, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> T, 
                                                   Eigen::Ref<const Eigen::MatrixXd> V_skin) {

    N.resize(3 * V_skin.rows(), 3 * V.rows());

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> triplets;

    for (int l = 0; l < V_skin.rows(); l++) {
        int temp_t = 0;
        double phi_norm = 999999999999999;
        Eigen::Vector4d phi;

        for (int i = 0; i < T.rows(); i++) {
            Eigen::Vector3d x = V_skin.row(l);
            Eigen::Vector4d temp_phi;
            Eigen::RowVectorXi element = T.row(i);
            phi_linear_tetrahedron(temp_phi, V, element, x);
            
            if (temp_phi.norm() < phi_norm) {
                temp_t = i;
                phi_norm = temp_phi.norm();
                phi = temp_phi;
            }
        }

        for(int z = 0; z < 4; z++){
            triplets.push_back(Triplet(l, T(temp_t, z), phi(z)));
        }
    }
    N.setFromTriplets(triplets.begin(), triplets.end());

}