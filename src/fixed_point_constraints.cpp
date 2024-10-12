#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

    std::vector<Eigen::Triplet<double> > temp_triplet;
    P.resize(q_size - 3 * indices.size(), q_size);
    P.setZero();

    int count = 0;
    for(int i = 0; i < q_size / 3; i++){
        if(count < indices.size() && indices[count] == i){count+=1;}else{
            for(int j = 0; j < 3; j++){
                temp_triplet.push_back(Eigen::Triplet<double>(3 * i - 3 * count + j, 3 * i + j, 1.0));
            }
        }
    }
    P.setFromTriplets(temp_triplet.begin(), temp_triplet.end());
}