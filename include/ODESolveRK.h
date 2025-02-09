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

    void alloc_func();
    void free_func();

    void init() {
        y_current_ = new double[n_func_y_ + 1];

        y_current_[0] = tab_func_y_[0][0];
        for (int i = 1; i <= n_func_y_; i++) {
            y_current_[i] = tab_func_y_[i][0];
        }

        initialized_ = true;

        return;
    }

    void evolve_RK4(double delta_x);
};

} // end namespace ODESolve

#endif