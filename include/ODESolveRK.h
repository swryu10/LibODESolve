#ifndef _ODESOLVERK_H_
#define _ODESOLVERK_H_

#include"ODESolveBase.h"

namespace ODESolve {

class LibRK : public LibBase {
  private :

    double *y_current_;

    bool have_func_rk_;

    bool initialized_;

  public :

    /* function pointers
     * to define differential equations
     *
     * For i, j = 1 ... n_func_y_,
     *   dy[i] / dx = g[i] (y[0 ... n_func_y_])
     * where y[0] = x. */
    func_ode_deriv *ptr_func_g_;

    LibRK() {
        ptr_log_ = NULL;

        have_func_rk_ = false;
        have_func_base_ = false;

        initialized_ = false;

        return;
    }

    ~LibRK() {
        if (!initialized_) {
            return;
        }

        delete [] y_current_;

        return;
    }

    /* The following variables must be specified
     * before calling alloc_func() function.
     *   n_func_y_
     *
     * After calling alloc_func() function,
     * one needs to specify the following variables.
     *   ptr_func_g_
     *
     * In addition, one has to feed the initial condition
     * into tab_func_y_[0 ... n_func_y_][0] */
    void alloc_func();

    /* a function to be called
     * at the end of program */
    void free_func();

    /* function to initialize the solver
     *
     * alloc_func() function must be called
     * and aforementioned variables have to be specified
     * before calling init() function. */
    void init() {
        y_current_ = new double[n_func_y_ + 1];

        y_current_[0] = tab_func_y_[0][0];
        for (int i = 1; i <= n_func_y_; i++) {
            y_current_[i] = tab_func_y_[i][0];
        }

        initialized_ = true;

        return;
    }

    /* function to implement 4-th order Runge-Kutta method
     *
     * Each time one calls evolve_RK4,
     * it appends y[i] at new x to tab_func_y_. */
    void evolve_RK4(double delta_x);
};

} // end namespace ODESolve

#endif