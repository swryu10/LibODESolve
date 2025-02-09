#include"ODESolveRK.h"

namespace ODESolve {

void LibRK::alloc_func() {
    alloc_func_base();
    for (int i = 0; i <= n_func_y_; i++) {
        tab_func_y_[i].resize(1);
    }

    ptr_func_g_ =
        (func_ode_deriv *)malloc((n_func_y_ + 1) *
                                 sizeof(func_ode_deriv));
    have_func_rk_ = true;

    return;
}

void LibRK::free_func() {
    if (!have_func_rk_) {
        return;
    }

    free_func_base();

    free(ptr_func_g_);

    return;
}

void LibRK::evolve_RK4(double delta_x) {
    if (!initialized_) {
        if (ptr_log_ != NULL) {
            fprintf(ptr_log_,
                "ODESolveRK ERROR : not initialized.\n");
        }

        exit(1);
    }

    double **y_comp = new double *[4];
    double **dy_dx_comp = new double *[4];
    for (int icomp = 0; icomp < 4; icomp++) {
        y_comp[icomp] = new double[n_func_y_ + 1];
        dy_dx_comp[icomp] = new double[n_func_y_ + 1];
    }

    y_comp[0][0] = y_current_[0];
    for (int i = 1; i <= n_func_y_; i++) {
        y_comp[0][i] = y_current_[i];
    }
    for (int i = 1; i <= n_func_y_; i++) {
        dy_dx_comp[0][i] =
            (*ptr_func_g_[i])(y_comp[0]);
    }

    y_comp[1][0] = y_current_[0] + 0.5 * delta_x;
    for (int i = 1; i <= n_func_y_; i++) {
        y_comp[1][i] =
            y_current_[i] + 0.5 * delta_x * dy_dx_comp[0][i];
    }
    for (int i = 1; i <= n_func_y_; i++) {
        dy_dx_comp[1][i] =
            (*ptr_func_g_[i])(y_comp[1]);
    }

    y_comp[2][0] = y_current_[0] + 0.5 * delta_x;
    for (int i = 1; i <= n_func_y_; i++) {
        y_comp[2][i] =
            y_current_[i] + 0.5 * delta_x * dy_dx_comp[1][i];
    }
    for (int i = 1; i <= n_func_y_; i++) {
        dy_dx_comp[2][i] =
            (*ptr_func_g_[i])(y_comp[2]);
    }

    y_comp[3][0] = y_current_[0] + delta_x;
    for (int i = 1; i <= n_func_y_; i++) {
        y_comp[3][i] =
            y_current_[i] + delta_x * dy_dx_comp[2][i];
    }
    for (int i = 1; i <= n_func_y_; i++) {
        dy_dx_comp[3][i] =
            (*ptr_func_g_[i])(y_comp[3]);
    }

    y_current_[0] += delta_x;
    tab_func_y_[0].push_back(y_current_[0]);
    for (int i = 1; i <= n_func_y_; i++) {
        y_current_[i] +=
            delta_x * (dy_dx_comp[0][i] +
                       2. * dy_dx_comp[1][i] +
                       2. * dy_dx_comp[2][i] +
                       dy_dx_comp[3][i]) / 6.;
        tab_func_y_[i].push_back(y_current_[i]);
    }

    for (int icomp = 0; icomp < 4; icomp++) {
        delete [] y_comp[icomp];
        delete [] dy_dx_comp[icomp];
    }
    delete [] y_comp;
    delete [] dy_dx_comp;

    return;
}

} // end namespace ODESolve
