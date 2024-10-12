 
 #include <mass_matrix_linear_tetrahedron.h>

 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {


     M.setZero();

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i == j){
                M.block(3 * i, 3 * j, 3, 3) = ((1.0 / 10.0) * density * volume) * Eigen::Matrix<double,3,3>::Identity();
            }else{
                M.block(3 * i, 3 * j, 3, 3) = ((1.0 / 20.0) * density * volume) * Eigen::Matrix<double,3,3>::Identity();
            }  
        }
    }
 }
