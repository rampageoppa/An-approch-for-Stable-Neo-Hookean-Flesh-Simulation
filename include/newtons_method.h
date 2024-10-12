#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  x0 - initial point for newtons search
//  f(x) - function that evaluates and returns the cost function at x
//  g(dx, x) - function that evaluates and returns the gradient of the cost function in dx
//  H(dH, x) - function that evaluates and returns the Hessian in dH (as a sparse matrix).
//  max steps - the maximum newton iterations to take
//  tmp_g and tmp_H are scratch space to store gradients and hessians
//Output: 
//  x0 - update x0 to new value
template<typename Objective, typename Jacobian, typename Hessian>
double newtons_method(Eigen::VectorXd &x0, Objective &f, Jacobian &g, Hessian &H, unsigned int maxSteps, Eigen::VectorXd &tmp_g, Eigen::SparseMatrixd &tmp_H) {

    for(int i = 0; i < maxSteps; i++)
    {
        g(tmp_g, x0);
        H(tmp_H, x0);
        if(tmp_g.norm() <= 1e-8)
        {
            break;
        }

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(tmp_H);
        Eigen::VectorXd d = solver.solve(-tmp_g);;

        double alpha = 1.0; 
        while(1){
            if(f(x0 + alpha * d) <= (f(x0) + 1e-8 * d.dot(tmp_g)) || alpha < 1e-8){    
                break;
            }
            alpha = 0.5 * alpha;
        }
        x0 = x0 + alpha * d;
    }

   return 0.0;
}
