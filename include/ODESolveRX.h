#ifndef _ODESOLVERX_H_
#define _ODESOLVERX_H_

#include"nrutil.h"
#include"ODESolveBase.h"

namespace ODESolve {

class LibRX : public LibBase {
  private :

    int ne_;
    int nb_;
    int mx_;

    long int nsi_;
    long int nsj_;
    long int nci_;
    long int ncj_;
    long int nck_;

    double *xbin_;
    double **mtx_y_;
    double **mtx_s_;
    double ***mtx_c_;

    int *kmax_;
    double *ermax_;

    int n_iteration_;

    bool have_func_rx_;

    bool initialized_;

  public :

    int n_boundary_i_;
    int n_bin_x_;

    double conv_;
    double slowc_;

    double *scalv_;
    int *indexv_;

    func_ode_deriv *ptr_func_g_;
    func_ode_deriv **ptr_func_dgdy_;

    func_ode_deriv *ptr_func_b_;
    func_ode_deriv **ptr_func_dbdy_;

    int n_iteration_max_;

    LibRX() {
        ptr_log_ = NULL;

        n_iteration_ = 0;
        n_iteration_max_ = 1000;

        have_func_rx_ = false;
        have_func_base_ = false;

        initialized_ = false;

        return;
    }

    ~LibRX() {
        free_array();

        return;
    }

    void alloc_func();
    void free_func();

    void free_array() {
        if (!initialized_) {
            return;
        }

        free(xbin_);
        free_dmatrix(mtx_y_, 1, ne_, 1, mx_);

        free_ivector(kmax_, 1, ne_);
        free_dvector(ermax_, 1, ne_);

        free_dmatrix(mtx_s_, 1, nsi_, 1, nsj_);
        free_d3tensor(mtx_c_, 1, nci_, 1, ncj_, 1, nck_);

        n_iteration_ = 0;

        initialized_ = false;

        return;
    }

    void init();
    bool next();

    void difeq_internal(int k, int k1, int k2,
                        int jsf, int is1, int isf);

    static bool solvde(int itmax, double conv, double slowc,
                       double *scalv, int *indexv,
                       int ne, int nb, int m,
                       double **y, double ***c, double **s,
                       void (*ptr_difeq)(int, int, int,
                                         int, int, int,
                                         int *, int,
                                         double **, double **));

    static void bksub(int ne, int nb, int jf,
                      int k1, int k2, double ***c);
    static void pinvs(int ie1, int ie2,
                      int je1, int jsf, int jc1,
                      int k, double ***c, double **s);
    static void red(int iz1, int iz2, int jz1, int jz2,
                    int jm1, int jm2, int jmf,
                    int ic1, int jc1, int jcf,
                    int kc, double ***c, double **s);
};

} // end namespace ODESolve

#endif