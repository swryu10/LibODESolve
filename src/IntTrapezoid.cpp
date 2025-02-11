#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#ifdef _OPENMP
#include<omp.h>
#endif
#include"IntTrapezoid.h"

namespace ODESolve {

FILE *ptr_log_IntTrapezoid_ = NULL;

double eps_precision_IntTrapezoid_ = 1e-8;

double get_integral_trapezoidal(double xmin, double xmax,
                                double (*ptr_func)(double)) {
    double ret_now = 0.;
    double ret_prev;
    int converging = 0;

    #ifdef _OPENMP
    #pragma omp parallel
    {  // parallel code begins
    #endif
        #ifdef _OPENMP
        int n_thread = omp_get_num_threads();
        int tid = omp_get_thread_num();
        if (ptr_log_IntTrapezoid_ != NULL) {
            fprintf(ptr_log_IntTrapezoid_,
                "  n_thread = %d, tid = %d\n",
                n_thread, tid);
        }
        #else
        int n_thread = 1;
        int tid = 0;
        #endif

        unsigned long int nbin_x = 1;
        double xmin_thread =
            xmin + (xmax - xmin) * (double)tid / (double)n_thread;
        double xmax_thread =
            xmin_thread + (xmax - xmin) / (double)n_thread;
        double delta_x =
            xmax_thread - xmin_thread;

        double ret_thread = 0.;
        double ret_thread_prev;

        int istep = 1;
        while (converging == 0) {
            if (tid == 0) {
                ret_prev = ret_now;
                ret_now = 0.;
            }

            #ifdef _OPENMP
            // syncronize threads
            #pragma omp barrier
            #endif

            ret_thread_prev = ret_thread;
            ret_thread = 0.;
            unsigned long int ix;
            if (istep == 1) {
                for (ix = 0; ix < nbin_x; ix++) {
                    double x0 = xmin_thread + delta_x * (double)ix;
                    double x1 = x0 + delta_x;
                    ret_thread +=
                        0.5 * delta_x * ((*ptr_func)(x0) +
                                         (*ptr_func)(x1));
                }
            } else {
                ret_thread = 0.5 * ret_thread_prev;
                for (ix = 0; ix <= nbin_x; ix++) {
                    if (ix % 2 == 0) {
                        continue;
                    }

                    double x_now = xmin_thread + delta_x * (double)ix;
                    ret_thread +=
                        delta_x * (*ptr_func)(x_now);
                }
            }

            #ifdef _OPENMP
            #pragma omp critical
            {  // exclusive block begins
            #endif
                ret_now += ret_thread;
            #ifdef _OPENMP
            }  // exclusive block ends
            #endif

            #ifdef _OPENMP
            // syncronize threads
            #pragma omp barrier
            #endif

            if (tid == 0) {
                if (ptr_log_IntTrapezoid_ != NULL) {
                    fprintf(ptr_log_IntTrapezoid_,
                        "      step %d : integral = %.10f\n",
                        istep, ret_now);
                }

                // conveergence check
                if (fabs(ret_now - ret_prev) <
                    0.5 * eps_precision_IntTrapezoid_ *
                    fabs(ret_now + ret_prev)) {
                    converging = 1;
                }
            }

            #ifdef _OPENMP
            // syncronize threads
            #pragma omp barrier
            #endif

            istep += 1;
            nbin_x = 2 * nbin_x;
            delta_x = 0.5 * delta_x;
        }
    #ifdef _OPENMP
    }  // parallel code ends
    #endif

    return ret_now;
}

} // end namespace ODESolve