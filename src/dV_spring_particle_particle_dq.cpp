#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {


    double norm = sqrt(pow(q0[0] - q1[0], 2.0) + pow(q0[1] - q1[1], 2.0) + pow(q0[2] - q1[2], 2.0));

    f[0] = stiffness * (l0 - norm) * (-1.0 / 2.0) * (q0[0] * 2.0 - q1[0] * 2.0) / norm;
    f[1] = stiffness * (l0 - norm) * (-1.0 / 2.0) * (q0[1] * 2.0 - q1[1] * 2.0) / norm;
    f[2] = stiffness * (l0 - norm) * (-1.0 / 2.0) * (q0[2] * 2.0 - q1[2] * 2.0) / norm;
    f[3] = stiffness * (l0 - norm) * (1.0 / 2.0) * (q0[0] * 2.0 - q1[0] * 2.0) / norm;
    f[4] = stiffness * (l0 - norm) * (1.0 / 2.0) * (q0[1] * 2.0 - q1[1] * 2.0) / norm;
    f[5] = stiffness * (l0 - norm) * (1.0 / 2.0) * (q0[2] * 2.0 - q1[2] * 2.0) / norm;

}