#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"IntTrapezoid.h"

double umin = 0.;
double umax = 1.;
double func_integrand(double u) {
    return 4. / (1. + u * u);
}

int main(int argc, char *argv[]) {
    ODESolve::ptr_log_IntTrapezoid_ = stderr;
    ODESolve::eps_precision_IntTrapezoid_ = 1e-8;

    double (*ptr_func_integrand)(double) = &func_integrand;

    double pi_now =
        ODESolve::get_integral_trapezoidal(umin, umax,
                                           ptr_func_integrand);

    //fprintf(stdout, "\n");
    fprintf(stdout, "pi from numerical integration\n");
    fprintf(stdout, "  > pi = %.8f\n", pi_now);
    fprintf(stdout, "pi from C math library\n");
    fprintf(stdout, "  > pi = %.8f\n", M_PI);

    return 0;
}
