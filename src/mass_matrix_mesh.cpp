#include <mass_matrix_mesh.h>
#include <mass_matrix_linear_tetrahedron.h>

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> T, double density, Eigen::Ref<const Eigen::VectorXd> v0) {

    M.setZero();
    M.resize(qdot.rows(), qdot.rows());

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> triplets;

    for(int i = 0; i < T.rows(); i++)
    {
        Eigen::Matrix1212d temp_m;
        mass_matrix_linear_tetrahedron(temp_m, qdot, T.row(i), density, v0(i));

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                for (int z = 0; z < 3; z++){
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k), temp_m(3 * j + z, 3 * k)));
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k) + 1, temp_m(3 * j + z, 3 * k + 1)));
                    triplets.push_back(Triplet(3 * T(i, j) + z, 3 * T(i, k) + 2, temp_m(3 * j + z, 3 * k + 2)));
                }
            }
        }     
    }
    M.setFromTriplets(triplets.begin(), triplets.end());
}
