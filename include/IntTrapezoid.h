#ifndef _INTTRAPEZOID_H_
#define _INTTRAPEZOID_H_

namespace ODESolve {

/* file pointer for
 * the log file, where we have
 * the approximate value at each step */
extern FILE *ptr_log_IntTrapezoid_;

/* precision to be acheived
 * If difference between the current step and previous one
 * is smaller than eps * (approximate) integration,
 * it is considered converging. */
extern double eps_precision_IntTrapezoid_;

/* prototype of the function,
 * which returns integral
 * from the trapezoidal rule
 * xmin : the lower limit of definite integral
 * xmax : the upper limit of definite integral
 * ptr_func : function pointer to the integrand */
double get_integral_trapezoidal(double xmin, double xmax,
                                double (*ptr_func)(double));

} // end namespace ODESolve

#endif
