#ifndef _ODESOLVEBASE_H_
#define _ODESOLVEBASE_H_

#include<stdio.h>
#include<stdlib.h>
#include<vector>

namespace ODESolve {

typedef double (*func_ode_deriv)(double *);

class LibBase {
  protected :

    bool have_func_base_;

  public :

    // file pointer to print event logs
    FILE *ptr_log_;

    // number of functions to be solved
    int n_func_y_;

    /* tabulated functions
     *
     * For i = 1 ... n_func_y_,
     * tab_func_y_[0][ix] = x at the ix-th mesh point
     * tab_func_y_[i][ix] = y[i] at the ix-th mesh point
     *
     * ix = 0 corresponds to the initial point
     * or the first boundary */
    std::vector<double> *tab_func_y_;

    LibBase() {
        ptr_log_ = NULL;

        have_func_base_ = false;

        return;
    }

    ~LibBase() {}

    void alloc_func_base();
    void free_func_base();

    void export_file(FILE *ptr_fout = stdout);
};

} // end namespacee ODESolve

#endif
