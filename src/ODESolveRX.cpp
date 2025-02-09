#include"ODESolveRX.h"

namespace ODESolve {

void LibRX::alloc_func() {
    alloc_func_base();
    for (int i = 0; i <= n_func_y_; i++) {
        tab_func_y_[i].resize(n_bin_x_ + 1);
    }

    scalv_ = new double[n_func_y_ + 1];
    indexv_ = new int[n_func_y_ + 1];

    ptr_func_g_ =
        (func_ode_deriv *)malloc((n_func_y_ + 1) *
                                 sizeof(func_ode_deriv));
    ptr_func_b_ =
        (func_ode_deriv *)malloc((n_func_y_ + 1) *
                                 sizeof(func_ode_deriv));

    ptr_func_dgdy_ =
        (func_ode_deriv **)malloc((n_func_y_ + 1) *
                                  sizeof(func_ode_deriv *));
    ptr_func_dbdy_ =
        (func_ode_deriv **)malloc((n_func_y_ + 1) *
                                  sizeof(func_ode_deriv *));

    for (int i = 1; i <= n_func_y_; i++) {
        ptr_func_dgdy_[i] =
            (func_ode_deriv *)malloc((n_func_y_ + 1) *
                                     sizeof(func_ode_deriv));
        ptr_func_dbdy_[i] =
            (func_ode_deriv *)malloc((n_func_y_ + 1) *
                                     sizeof(func_ode_deriv));
    }

    have_func_rx_ = true;

    return;
}

void LibRX::free_func() {
    if (!have_func_rx_) {
        return;
    }

    free_func_base();

    delete [] scalv_;
    delete [] indexv_;

    for (int j = 0; j < n_func_y_; j++) {
        free(ptr_func_dgdy_[j]);
        free(ptr_func_dbdy_[j]);
    }
    free(ptr_func_dgdy_);
    free(ptr_func_dbdy_);

    free(ptr_func_g_);
    free(ptr_func_b_);

    return;
}

void LibRX::init() {
    free_array();

    ne_ = n_func_y_;
    nb_ = n_boundary_i_;
    mx_ = n_bin_x_ + 1;

    xbin_ = (double *)malloc((mx_ + 1) * sizeof(double));
    mtx_y_ = dmatrix(1, ne_, 1, mx_);
    for (int ix = 1; ix <= mx_; ix++) {
        xbin_[ix] = tab_func_y_[0][ix - 1];

        for (int i = 1; i <= ne_; i++) {
            mtx_y_[i][ix] = tab_func_y_[i][ix - 1];
        }
    }

    kmax_ = ivector(1, ne_);
    ermax_ = dvector(1, ne_);

    nsi_ = ne_;
    nsj_ = 2 * ne_ + 1;
    mtx_s_ = dmatrix(1, nsi_, 1, nsj_);

    nci_ = ne_;
    ncj_ = ne_ - nb_ + 1;
    nck_ = mx_ + 1;
    mtx_c_ = d3tensor(1, nci_, 1, ncj_, 1, nck_);

    n_iteration_ = 0;

    initialized_ = true;

    return;
}

bool LibRX::next() {
    if (!initialized_) {
        if (ptr_log_ != NULL) {
            fprintf(ptr_log_,
                "ODESolveRX ERROR : not initialized.\n");
        }

        exit(1);
    }

    if (n_iteration_ >= n_iteration_max_) {
        if (ptr_log_ != NULL) {
            fprintf(ptr_log_,
                "ODESolveRX WARNING : nmax_iteration reached.\n");
        }

        return true;
    }

    int k1 = 1;
    int k2 = mx_;
    int nvars = ne_ * mx_;
    int j1 = 1;
    int j2 = nb_;
    int j3 = nb_ + 1;
    int j4 = ne_;
    int j5 = j4 + j1;
    int j6 = j4 + j2;
    int j7 = j4 + j3;
    int j8 = j4 + j4;
    int j9 = j8 + j1;
    int ic1 = 1;
    int ic2 = ne_ - nb_;
    int ic3 = ic2 + 1;
    int ic4 = ne_;
    int jc1 = 1;
    int jcf = ic3;

    //Boundary conditions at first point.
    int k = k1;
    difeq_internal(k, k1, k2, j9, ic3, ic4);
    pinvs(ic3, ic4, j5, j9, jc1, k1, mtx_c_, mtx_s_);

    for (k = k1 + 1; k <= k2; k++) {
        /* Finite difference equations
         * at all point pairs. */
        int kp = k - 1;
        difeq_internal(k, k1, k2, j9, ic1, ic4);
        red(ic1, ic4, j1, j2, j3, j4, j9,
            ic3, jc1, jcf, kp, mtx_c_, mtx_s_);
        pinvs(ic1, ic4, j3, j9, jc1, k, mtx_c_, mtx_s_);
    }

    k = k2 + 1;  //Final boundary conditions.
    difeq_internal(k, k1, k2, j9, ic1, ic2);
    red(ic1, ic2, j5, j6, j7, j8, j9,
        ic3, jc1, jcf, k2, mtx_c_, mtx_s_);
    pinvs(ic1, ic2, j7, j9, jcf, k2 + 1, mtx_c_, mtx_s_);

    //Backsubstitution.
    bksub(ne_, nb_, jcf, k1, k2, mtx_c_);
    double err = 0.;

    for (int j = 1; j <= ne_; j++) {
        /* Convergence check
         * accumulate average error */
        int jv = indexv_[j];
        int km = 0;
        double vmax = 0.;
        double errj = vmax;
        for (k = k1; k <= k2; k++) {
            /* Find point with largest error,
             * for each dependent variable. */
            double vz = fabs(mtx_c_[jv][1][k]);
            if (vz > vmax) {
                vmax = vz;
                km = k;
            }
            errj += vz;
        }
        err += errj / scalv_[j];
        // Note weighting for each dependent variable.
        ermax_[j] = mtx_c_[jv][1][km] / scalv_[j];

        kmax_[j] = km;
    }

    err /= (double)nvars;
    double fac =
        (err > slowc_ ? slowc_ / err : 1.0);
    //Reduce correction applied when error is large.
    for (int j = 1; j <= ne_; j++) {  // Apply corrections.
        int jv = indexv_[j];
        for (k = k1; k <= k2; k++) {
            mtx_y_[j][k] -= fac * mtx_c_[jv][1][k];
        }
    }

    n_iteration_ += 1;

    if (ptr_log_ != NULL) {
        fprintf(ptr_log_, "\n");
        /* Summary of corrections
         * for this step */
        fprintf(ptr_log_,
            "    n_iteration = %d\n", n_iteration_);
        fprintf(ptr_log_,
            "      err = %12.6f, fac = %11.6f\n", err, fac);
    }

    bool have_final = false;

    if (err < conv_) {
        have_final = true;

        for (int ix = 1; ix <= mx_; ix++) {
            for (int i = 1; i <= ne_; i++) {
                tab_func_y_[i][ix - 1] = mtx_y_[i][ix];
            }
        }
    }

    return have_final;
}

void LibRX::difeq_internal(int k, int k1, int k2,
                           int jsf, int is1, int isf) {
    double *ymed;
    ymed = (double *)malloc((ne_ + 1) * sizeof(double));

    //fprintf(stderr, "\n");
    //fprintf(stderr, "      ymed[] = ");
    for (int i = 1; i <= ne_; i++) {
        if (k == k1) {
            ymed[i] = mtx_y_[i][k1];
        } else if (k > k2) {
            ymed[i] = mtx_y_[i][k2];
        } else {
            ymed[i] = 0.5 * (mtx_y_[i][k] + mtx_y_[i][k - 1]);
        }

        //fprintf(stderr, "  %e", ymed[i]);
    }
    //fprintf(stderr, "\n");

    if (k == k1) {
        // Boundary condition at the first point
        ymed[0] = xbin_[k1];

        for (int i = ne_ - nb_ + 1; i <= ne_; i++) {
            if (i < is1 || i > isf) {
                continue;
            }

            for (int j = 1; j <= ne_; j++) {
                int jv = indexv_[j];

                mtx_s_[i][jv] = 0.;
                mtx_s_[i][jv + ne_] =
                    (*ptr_func_dbdy_[i - ne_ + nb_][j])(ymed);
            }

            mtx_s_[i][jsf] =
                (*ptr_func_b_[i - ne_ + nb_])(ymed);
        }
    } else if (k > k2) {
        // Boundary condition at the last point
        ymed[0] = xbin_[k2];

        for (int i = 1; i <= ne_ - nb_; i++) {
            if (i < is1 || i > isf) {
                continue;
            }

            for (int j = 1; j <= ne_; j++) {
                int jv = indexv_[j];

                mtx_s_[i][jv] = 0.;
                mtx_s_[i][jv + ne_] =
                    (*ptr_func_dbdy_[i + nb_][j])(ymed);
            }

            mtx_s_[i][jsf] =
                (*ptr_func_b_[i + nb_])(ymed);
        }
    } else {
        // Interior points
        ymed[0] = 0.5 * (xbin_[k] + xbin_[k - 1]);

        double xdiff = xbin_[k] - xbin_[k - 1];

        for (int i = 1; i <= ne_; i++) {
            if (i < is1 || i > isf) {
                continue;
            }

            for (int j = 1; j <= ne_; j++) {
                int jv = indexv_[j];

                mtx_s_[i][jv] =
                    -0.5 * xdiff *
                    (*ptr_func_dgdy_[i][j])(ymed);
                mtx_s_[i][jv + ne_] = mtx_s_[i][jv];
                if (i == j) {
                    mtx_s_[i][jv] -= 1.;
                    mtx_s_[i][jv + ne_] += 1.;
                }
            }

            mtx_s_[i][jsf] =
                mtx_y_[i][k] - mtx_y_[i][k - 1] -
                xdiff * (*ptr_func_g_[i])(ymed);
        }
    }

    /*
    fprintf(stderr, "\n");
    for (int i = 1; i <= ne_; i++) {
        if (i < is1 || i > isf) {
            continue;
        }

        fprintf(stderr, "    i = %d :", i);
        for (int j = 1; j <= 2 * ne_ + 1; j++) {
            if (j < 1 || j > jsf) {
                continue;
            }

            fprintf(stderr, "  %e", mtx_s_[i][j]);
        }
        fprintf(stderr, "\n");
    }
    */

    free(ymed);
}

bool LibRX::solvde(int itmax, double conv, double slowc,
                   double *scalv, int *indexv,
                   int ne, int nb, int m,
                   double **y, double ***c, double **s,
    void (*ptr_difeq)(int, int, int, int, int, int,
                      int *, int,
                      double **, double **)) {
    /* Driver routine for solution
     * of two point boundary value problem by relaxation.
     * itmax is the maximum number of iterations.
     * conv is the convergence criterion (see text).
     * slowc controls the fraction of corrections actually used
     * after each iteration.
     * scalv[1..ne] contains typical sizes
     * for each dependent variable, used to weight errors.
     * indexv[1..ne] lists the column ordering of variables
     * used to construct the matrix
     * s[1..ne][1..2*ne+1] of derivatives.
     * (The nb boundary conditions at the first mesh point
     *  must contain some dependence on the first nb
     *  variables listed in indexv.)
     * The problem involves ne equations
     * for ne adjustable dependent variables at each point.
     * At the first mesh point there are nb boundary conditions.
     * There are a total of m mesh points.
     * y[1..ne][1..m] is the two-dimensional array
     * that contains the initial guess
     * for all the dependent variables at each mesh point.
     * On each iteration, it is updated by the calculated correction.
     * The arrays c[1..ne][1..ne-nb+1][1..m+1] and s supply
     * dummy storage used by the relaxation code. */

    int *kmax = ivector(1, ne);
    double *ermax = dvector(1, ne);

    int k1 = 1;
    int k2 = m;
    int nvars = ne * m;
    int j1 = 1;
    int j2 = nb;
    int j3 = nb + 1;
    int j4 = ne;
    int j5 = j4 + j1;
    int j6 = j4 + j2;
    int j7 = j4 + j3;
    int j8 = j4 + j4;
    int j9 = j8 + j1;
    int ic1 = 1;
    int ic2 = ne - nb;
    int ic3 = ic2 + 1;
    int ic4 = ne;
    int jc1 = 1;
    int jcf = ic3;

    bool converging = false;

    // Primary iteration loop.
    for (int it = 1; it <= itmax; it++) {
        //Boundary conditions at first point.
        int k = k1;
        (*ptr_difeq)(k, k1, k2, j9, ic3, ic4,
                     indexv, ne, s, y);
        pinvs(ic3, ic4, j5, j9, jc1, k1, c, s);
        for (k = k1 + 1; k <= k2; k++) {
            /* Finite difference equations
             * at all point pairs. */
            int kp = k - 1;
            (*ptr_difeq)(k, k1, k2, j9, ic1, ic4,
                         indexv, ne, s, y);
            red(ic1, ic4, j1, j2, j3, j4, j9,
                ic3, jc1, jcf, kp, c, s);
            pinvs(ic1, ic4, j3, j9, jc1, k, c, s);
        }
        k = k2 + 1;  //Final boundary conditions.
        (*ptr_difeq)(k, k1, k2, j9, ic1, ic2,
                     indexv, ne, s, y);
        red(ic1, ic2, j5, j6, j7, j8, j9,
            ic3, jc1, jcf, k2, c, s);
        pinvs(ic1, ic2, j7, j9, jcf, k2 + 1, c, s);

        //Backsubstitution.
        bksub(ne, nb, jcf, k1, k2, c);
        double err = 0.;

        for (int j = 1; j <= ne; j++) {
            /* Convergence check
             * accumulate average error */
            int jv = indexv[j];
            int km = 0;
            double vmax = 0.;
            double errj = vmax;
            for (k = k1; k <= k2; k++) {
                /* Find point with largest error,
                 * for each dependent variable. */
                double vz = fabs(c[jv][1][k]);
                if (vz > vmax) {
                    vmax = vz;
                    km = k;
                }
                errj += vz;
            }
            err += errj / scalv[j];
            // Note weighting for each dependent variable.
            ermax[j] = c[jv][1][km] / scalv[j];

            kmax[j] = km;
        }

        err /= (double)nvars;
        double fac =
            (err > slowc ? slowc / err : 1.0);
        //Reduce correction applied when error is large.
        for (int j = 1; j <= ne; j++) {  // Apply corrections.
            int jv = indexv[j];
            for (k = k1; k <= k2; k++) {
                y[j][k] -= fac * c[jv][1][k];
            }
        }
        /* Summary of corrections
         * for this step */
        printf("\n%8s %9s %9s\n", "Iter.", "Error", "FAC");
        printf("%6d %12.6f %11.6f\n", it, err, fac);
        if (err < conv) {
            /* Point with largest error for each variable
             * can be monitored
             * by writing out kmax_ and ermax_. */

            converging = true;
            break;
        }
    }

    free_ivector(kmax, 1, ne);
    free_dvector(ermax, 1, ne);

    return converging;
}

void LibRX::bksub(int ne, int nb, int jf,
                  int k1, int k2, double ***c) {
    /* Backsubstitution,
     * used internally by solvde. */
    int nbf = ne - nb;
    int im = 1;

    for (int k = k2; k >= k1; k--) {
        /* Use recurrence relations to eliminate
         * remaining dependences. */
        if (k == k1) {
            im = nbf + 1;
        }

        int kp = k + 1;

        for (int j = 1; j <= nbf; j++) {
            double xx = c[j][jf][kp];
            for (int i = im; i <= ne; i++) {
                c[i][jf][k] -= c[i][j][k] * xx;
            }
        }
    }

    for (int k = k1; k <= k2; k++) {
        // Reorder corrections to be in column 1.
        int kp = k + 1;

        for (int i = 1; i <= nb; i++) {
            c[i][1][k] = c[i + nbf][jf][k];
        }

        for (int i = 1; i <= nbf; i++) {
            c[i + nb][1][k] = c[i][jf][kp];
        }
    }
}

void LibRX::pinvs(int ie1, int ie2,
                  int je1, int jsf, int jc1,
                  int k, double ***c, double **s) {
    /* Diagonalize the square subsection
     * of the s matrix,
     * and store the recursion coefficients in c
     * used internally by solvde. */

    int jpiv, jp, jcoff, j, irow, ipiv, id, icoff, i;
    double pivinv, piv, dum, big;

    int *indxr = ivector(ie1,ie2);
    double *pscl = dvector(ie1,ie2);

    int je2 = je1 + ie2 - ie1;
    int js1 = je2 + 1;

    for (i = ie1; i <= ie2; i++) {
        // Implicit pivoting, as in section 2.1.
        big = 0.;
        for (j = je1; j <= je2; j++) {
            if (fabs(s[i][j]) > big) {
                big = fabs(s[i][j]);
            }
        }
        if (big == 0.) {
            nrerror("Singular matrix - row all 0, in pinvs");
        }
        pscl[i] = 1. / big;
        indxr[i] = 0;
    }

    for (id = ie1; id <= ie2; id++) {
        piv = 0.;

        for (i = ie1; i <= ie2; i++) {
            // Find pivot element.
            if (indxr[i] == 0) {
                big = 0.;
                for (j = je1; j <= je2; j++) {
                    if (fabs(s[i][j]) > big) {
                        jp = j;
                        big = fabs(s[i][j]);
                    }
                }
                if (big * pscl[i] > piv) {
                    ipiv = i;
                    jpiv = jp;
                    piv = big * pscl[i];
                }
            }
        }

        if (s[ipiv][jpiv] == 0.) {
            nrerror("Singular matrix in routine pinvs");
        }

        // In place reduction. Save column ordering.
        indxr[ipiv] = jpiv;
        pivinv = 1. / s[ipiv][jpiv];
        for (j = je1; j <= jsf; j++) {
            // Normalize pivot row.
            s[ipiv][j] *= pivinv;
        }
        s[ipiv][jpiv] = 1.;
        for (i = ie1; i <= ie2; i++) {
            //Reduce nonpivot elements in column.
            if (indxr[i] != jpiv) {
                if (s[i][jpiv]) {
                    dum = s[i][jpiv];
                    for (j = je1; j <= jsf; j++) {
                        s[i][j] -= dum * s[ipiv][j];
                    }
                    s[i][jpiv] = 0.;
                }
            }
        }
    }

    // Sort and store unreduced coefficients.
    jcoff = jc1 - js1;
    icoff = ie1 - je1;
    for (i = ie1; i <= ie2; i++) {
        irow = indxr[i] + icoff;
        for (j = js1; j <= jsf; j++) {
            c[irow][j + jcoff][k] = s[i][j];
        }
    }

    free_ivector(indxr, ie1, ie2);
    free_dvector(pscl, ie1, ie2);
}

void LibRX::red(int iz1, int iz2, int jz1, int jz2,
                int jm1, int jm2, int jmf,
                int ic1, int jc1, int jcf,
                int kc, double ***c, double **s) {
/* Reduce columns jz1 - jz2 of the s matrix,
 * using previous results as stored in the c matrix.
 * Only columns jm1 - jm2, jmf are affected by the prior results.
 * red is used internally by solvde. */
    int l, j, i;
    double vx;
    int loff = jc1 - jm1;
    int ic = ic1;

    for (j = jz1; j <= jz2; j++) {
        // Loop over columns to be zeroed.
        for (l = jm1; l <= jm2; l++) {
            // Loop over columns altered.
            vx = c[ic][l + loff][kc];
            for (i = iz1; i <= iz2; i++) {
                // Loop over rows.
                s[i][l] -= s[i][j] * vx;
            }
        }
        vx = c[ic][jcf][kc];
        for (i = iz1; i <= iz2; i++) {
            // Plus final element.
            s[i][jmf] -= s[i][j] * vx;
        }
        ic += 1;
    }
}

} // end namespace ODESolve
