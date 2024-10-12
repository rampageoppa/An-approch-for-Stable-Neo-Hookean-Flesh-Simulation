#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    // V = 1/2k(l-l_0)^2
    V = 0.5 * stiffness * pow((q0 - q1).norm() - l0, 2);
    
}