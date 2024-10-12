#include <dpsi_neo_hookean_dF.h>

void dpsi_neo_hookean_dF(Eigen::Vector9d &dw, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D, bool stable) {

    double F1_1 = F(0, 0);
    double F1_2 = F(0, 1);
    double F1_3 = F(0, 2);
    double F2_1 = F(1, 0);
    double F2_2 = F(1, 1);
    double F2_3 = F(1, 2);
    double F3_1 = F(2, 0);
    double F3_2 = F(2, 1);
    double F3_3 = F(2, 2);

    if(stable)
    {
        dw[0] = F1_1 * (2.0 * C) - (2.0 * D) * (F2_2 * F3_3 - F2_3 * F3_2) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[1] = F1_2 * (2.0 * C) + (2.0 * D) * (F2_1 * F3_3 - F2_3 * F3_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[2] = F1_3 * (2.0 * C) - (2.0 * D) * (F2_1 * F3_2 - F2_2 * F3_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[3] = F2_1 * (2.0 * C) + (2.0 * D) * (F1_2 * F3_3 - F1_3 * F3_2) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[4] = F2_2 * (2.0 * C) - (2.0 * D) * (F1_1 * F3_3 - F1_3 * F3_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[5] = F2_3 * (2.0 * C) + (2.0 * D) * (F1_1 * F3_2 - F1_2 * F3_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[6] = F3_1 * (2.0 * C) - (2.0 * D) * (F1_2 * F2_3 - F1_3 * F2_2) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[7] = F3_2 * (2.0 * C) + (2.0 * D) * (F1_1 * F2_3 - F1_3 * F2_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
        dw[8] = F3_3 * (2.0 * C) - (2.0 * D) * (F1_1 * F2_2 - F1_2 * F2_1) * ((2.0 * C) / (2.0 * D) - F1_1 * F2_2 * F3_3 + F1_1 * F2_3 * F3_2 + F1_2 * F2_1 * F3_3 - F1_2 * F2_3 * F3_1 - F1_3 * F2_1 * F3_2 + F1_3 * F2_2 * F3_1 + 1.0);
    }
    else{
        dw(0,0) = C * (F1_1 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 - (F2_2 * F3_3 - F2_3 * F3_2) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) - D * (F2_2 * F3_3 - F2_3 * F3_2) * ( - F.determinant() + 1.0) * 2.0;
        dw(1,0) = C * (F1_2 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 + (F2_1 * F3_3 - F2_3 * F3_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) + D * (F2_1 * F3_3 - F2_3 * F3_1) * ( - F.determinant() + 1.0) * 2.0;
        dw(2,0) = C * (F1_3 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 - (F2_1 * F3_2 - F2_2 * F3_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) - D * (F2_1 * F3_2 - F2_2 * F3_1) * ( - F.determinant() + 1.0) * 2.0;
        dw(3,0) = C * (F2_1 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 + (F1_2 * F3_3 - F1_3 * F3_2) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) + D * (F1_2 * F3_3 - F1_3 * F3_2) * ( - F.determinant() + 1.0) * 2.0;
        dw(4,0) = C * (F2_2 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 - (F1_1 * F3_3 - F1_3 * F3_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) - D * (F1_1 * F3_3 - F1_3 * F3_1) * ( - F.determinant() + 1.0) * 2.0;
        dw(5,0) = C * (F2_3 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 + (F1_1 * F3_2 - F1_2 * F3_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) + D * (F1_1 * F3_2 - F1_2 * F3_1) * ( - F.determinant() + 1.0) * 2.0;
        dw(6,0) = C * (F3_1 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 - (F1_2 * F2_3 - F1_3 * F2_2) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) - D * (F1_2 * F2_3 - F1_3 * F2_2) * ( - F.determinant() + 1.0) * 2.0;
        dw(7,0) = C * (F3_2 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 + (F1_1 * F2_3 - F1_3 * F2_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) + D * (F1_1 * F2_3 - F1_3 * F2_1) * ( - F.determinant() + 1.0) * 2.0;
        dw(8,0) = C * (F3_3 * 1.0 / pow(F.determinant(), 2.0 / 3.0) * 2.0 - (F1_1 * F2_2 - F1_2 * F2_1) * 1.0 / pow(F.determinant(), 5.0 / 3.0) * (F.transpose() * F).trace() * (2.0 / 3.0)) - D * (F1_1 * F2_2 - F1_2 * F2_1) * ( - F.determinant() + 1.0) * 2.0;
    }
}