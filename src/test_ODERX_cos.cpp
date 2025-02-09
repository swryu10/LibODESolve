#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"ODESolveRX.h"

double y2_ini = 0.;
double y1_fin = 1.;

double func_g1(double *y) {
    return y[2];
}
double func_g2(double *y) {
    return -y[1];
}

double func_b1(double *y) {
    return y[2] - y2_ini;
}
double func_b2(double *y) {
    return y[1] - y1_fin;
}

double func_dg11(double *y) {
    return 0.;
}
double func_dg12(double *y) {
    return 1.;
}
double func_dg21(double *y) {
    return -1.;
}
double func_dg22(double *y) {
    return 0.;
}

double func_db11(double *y) {
    return 0.;
}
double func_db12(double *y) {
    return 1.;
}
double func_db21(double *y) {
    return 1.;
}
double func_db22(double *y) {
    return 0.;
}

int main(int argc, char *argv[]) {
    ODESolve::LibRX odesol_rx;

    odesol_rx.ptr_log_ = stderr;

    odesol_rx.n_func_y_ = 2;
    odesol_rx.n_boundary_i_ = 1;

    double xmin = 0.;
    double xmax = 2. * M_PI;
    double delta_x = 0.02;
    odesol_rx.n_bin_x_ =
        (int)ceil(fabs(xmax - xmin) / delta_x);
    delta_x = fabs(xmax - xmin) /
              static_cast<double>(odesol_rx.n_bin_x_);

    odesol_rx.alloc_func();
    for (int ix = 0; ix <= odesol_rx.n_bin_x_; ix++) {
        odesol_rx.tab_func_y_[0][ix] = xmin + delta_x * (double)ix;

        odesol_rx.tab_func_y_[1][ix] = 1.;
        odesol_rx.tab_func_y_[2][ix] = 0.;
    }

    odesol_rx.conv_ = 1.0e-6;
    odesol_rx.slowc_ = 0.8;

    odesol_rx.scalv_[1] = 1.;
    odesol_rx.scalv_[2] = 1.;

    odesol_rx.indexv_[1] = 2;
    odesol_rx.indexv_[2] = 1;

    odesol_rx.ptr_func_g_[1] = &func_g1;
    odesol_rx.ptr_func_g_[2] = &func_g2;

    odesol_rx.ptr_func_b_[1] = &func_b1;
    odesol_rx.ptr_func_b_[2] = &func_b2;

    odesol_rx.ptr_func_dgdy_[1][1] = &func_dg11;
    odesol_rx.ptr_func_dgdy_[1][2] = &func_dg12;
    odesol_rx.ptr_func_dgdy_[2][1] = &func_dg21;
    odesol_rx.ptr_func_dgdy_[2][2] = &func_dg22;

    odesol_rx.ptr_func_dbdy_[1][1] = &func_db11;
    odesol_rx.ptr_func_dbdy_[1][2] = &func_db12;
    odesol_rx.ptr_func_dbdy_[2][1] = &func_db21;
    odesol_rx.ptr_func_dbdy_[2][2] = &func_db22;

    odesol_rx.init();

    bool converging = false;
    while (!converging) {
        converging = odesol_rx.next();
    }

    FILE *fout;
    fout = fopen("tab_oderx_cos.txt", "w");

    odesol_rx.export_file(fout);
    odesol_rx.free_func();

    fclose(fout);

    return 0;
}
