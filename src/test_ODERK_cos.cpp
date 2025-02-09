#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"ODESolveRK.h"

double y1_ini = 1.;
double y2_ini = 0.;

double func_g1(double *y) {
    return y[2];
}
double func_g2(double *y) {
    return -y[1];
}

int main(int argc, char *argv[]) {
    ODESolve::LibRK odesol_rk;

    odesol_rk.ptr_log_ = stderr;

    odesol_rk.n_func_y_ = 2;

    double xmin = 0.;
    double xmax = 2. * M_PI;
    double delta_x = 0.02;
    int n_bin_x =
        (int)ceil(fabs(xmax - xmin) / delta_x);
    delta_x = fabs(xmax - xmin) /
              static_cast<double>(n_bin_x);

    odesol_rk.alloc_func();
    odesol_rk.tab_func_y_[0][0] = xmin;
    odesol_rk.tab_func_y_[1][0] = y1_ini;
    odesol_rk.tab_func_y_[2][0] = y2_ini;

    odesol_rk.ptr_func_g_[1] = &func_g1;
    odesol_rk.ptr_func_g_[2] = &func_g2;

    odesol_rk.init();

    for (int ix = 0; ix < n_bin_x; ix++) {
        odesol_rk.evolve_RK4(delta_x);
    }

    FILE *fout;
    fout = fopen("tab_oderk_cos.txt", "w");

    odesol_rk.export_file(fout);
    odesol_rk.free_func();

    fclose(fout);

    return 0;
}
