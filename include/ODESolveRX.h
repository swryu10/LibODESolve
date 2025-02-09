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

    /* number of boundary conditions
     * at the first boundary
     *
     * There are n_func_y_ - n_boundary_i_ boundary conditions
     * at the second boundary. */
    int n_boundary_i_;

    /* number of bins in the x-direction
     * number of mesh points becomes n_bin_x_ + 1 */
    int n_bin_x_;

    /* parameter
     * to set precision for convergence check */
    double conv_;

    /* parameter
     * to modulate application of corrections */
    double slowc_;

    /* typical size
     * of each function */
    double *scalv_;

    /* specify ordering of functions
     * such that there is no zero pivot element in s[i][j]
     *
     * For instance, if there are three functions and
     * y[2] and y[3] are specified at the first boundary,
     * we have
     *   indexv_[1] = 3
     *   indexv_[2] = 1
     *   indexv_[3] = 2 */
    int *indexv_;

    /* function pointers
     * to define differential equations
     *
     * For i, j = 1 ... n_func_y_,
     *   dy[i] / dx = g[i] (y[0 ... n_func_y_])
     *   dgdy[i][j] = dg[i] / dy[j] (y[0 ... n_func_y_])
     * where y[0] = x. */
    func_ode_deriv *ptr_func_g_;
    func_ode_deriv **ptr_func_dgdy_;

    /* function pointers
     * to define boundary conditions
     *
     * For i, j = 1 ... n_func_y_,
     *   b[i] (y[0 ... n_func_y_]) = 0
     *   dbdy[i][j] = db[i] / dy[j] (y[0 ... n_func_y_])
     * where y[0] = x.
     *
     * i = 1 ... n_boundary_i_
     * are reserved for the conditions at the first boundary.
     * i = n_boundary_i_ + 1 ... n_func_y_
     * are reserved for the conditions at the second boundary. */
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

    /* The following variables must be specified
     * before calling alloc_func() function.
     *   n_func_y_
     *   n_boundary_i_
     *   n_bin_x_
     *
     * After calling alloc_func() function,
     * one needs to specify the following variables.
     *   conv_
     *   slowc_
     *   scalv_
     *   indexv_
     *   ptr_func_g_
     *   ptr_func_b_
     *   ptr_func_dgdy_
     *   ptr_func_dbdy_
     *
     * In addition, one has to feed the initial guess
     * into tab_func_y_[0 ... n_func_y_][0 ... n_bin_x_] */
    void alloc_func();

    /* a function to be called
     * at the end of program */
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

    /* function to initialize the solver
     *
     * alloc_func() function must be called
     * and aforementioned variables have to be specified
     * before calling init() function. */
    void init();

    /* function to implement relaxation method
     *
     * One needs to repeatedly call next() function
     * until it returns true. */
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