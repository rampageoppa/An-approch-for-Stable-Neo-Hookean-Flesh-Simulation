#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {

    //dense types
    using Vector4d = Eigen::Matrix<double,4,1>;
    using Vector6d = Eigen::Matrix<double, 6,1>;
    using Vector9d = Eigen::Matrix<double, 9, 1>;
    using Vector12d = Eigen::Matrix<double, 12,1>;
    using Matrix36d = Eigen::Matrix<double, 3,6>;
    using Matrix34d = Eigen::Matrix<double, 3,4>;
    using Matrix43d = Eigen::Matrix<double, 4,3>;
    using Matrix66d  = Eigen::Matrix<double, 6,6>;
    using Matrix99d = Eigen::Matrix<double, 9, 9>;
    using Matrix1212d = Eigen::Matrix<double, 12,12>;
    using Matrix44f = Eigen::Matrix<float, 4,4>;
    
    //sparse types
    using SparseMatrixd = Eigen::SparseMatrix<double>;

}

#endif 
