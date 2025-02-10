#include"OrbitKepler.h"

double l_ang = 1.;

double func_g1(double *y) {
    return y[2];
}
double func_g2(double *y) {
    return l_ang * l_ang / (y[1] * y[1] * y[1]) -
           1. / (y[1] * y[1]);
}
double func_g3(double *y) {
    return l_ang / (y[1] * y[1]);
}

void OrbitKepler::init(double r_i_in,
                       double vr_i_in,
                       double phi_i_in,
                       double l_angular_in,
                       double t_min_in,
                       double t_max_in,
                       double delta_t_in) {
    free_odesolve();

    r_ini_ = r_i_in;
    vr_ini_ = vr_i_in;
    phi_ini_ = phi_i_in;

    l_angular_ = l_angular_in;
    l_ang = l_angular_in;

    ptr_odesol_ = new ODESolve::LibRK();

    ptr_odesol_->ptr_log_ = stderr;

    ptr_odesol_->n_func_y_ = 3;

    t_min_ = t_min_in;
    t_max_ = t_max_in;
    delta_t_ = delta_t_in;
    n_bin_t_ =
        (int)ceil(fabs(t_max_ - t_min_) / delta_t_);
    delta_t_ = fabs(t_max_ - t_min_) /
              static_cast<double>(n_bin_t_);

    ptr_odesol_->alloc_func();
    ptr_odesol_->tab_func_y_[0][0] = t_min_;
    ptr_odesol_->tab_func_y_[1][0] = r_ini_;
    ptr_odesol_->tab_func_y_[2][0] = vr_ini_;
    ptr_odesol_->tab_func_y_[3][0] = phi_ini_;

    ptr_odesol_->ptr_func_g_[1] = &func_g1;
    ptr_odesol_->ptr_func_g_[2] = &func_g2;
    ptr_odesol_->ptr_func_g_[3] = &func_g3;

    ptr_odesol_->init();

    tab_t_.resize(n_bin_t_ + 1);
    tab_x_.resize(n_bin_t_ + 1);
    tab_y_.resize(n_bin_t_ + 1);

    initialized_ = true;

    return;
}

void OrbitKepler::find_orbit() {
    if (!initialized_) {
        return;
    }

    for (int it = 0; it < n_bin_t_; it++) {
        ptr_odesol_->evolve_RK4(delta_t_);
    }

    for (int it = 0; it <= n_bin_t_; it++) {
        tab_t_[it] = ptr_odesol_->tab_func_y_[0][it];

        double r_now = ptr_odesol_->tab_func_y_[1][it];
        double phi_now = ptr_odesol_->tab_func_y_[3][it];

        tab_x_[it] = r_now * cos(phi_now);
        tab_y_[it] = r_now * sin(phi_now);
    }

    return;
}
