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

    FILE *ptr_log_;

    int n_func_y_;
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
