#include"GeodesicSphere.h"

double phi_rad_i = 0.;
double phi_rad_f = 0.;

double lambda_rad_i = 0.;
double lambda_rad_f = 0.;

double func_g1(double *y) {
    return y[3];
}
double func_g2(double *y) {
    return y[4];
}
double func_g3(double *y) {
    return -cos(y[1]) * sin(y[1]) * y[4] * y[4];
}
double func_g4(double *y) {
    return 2. * tan(y[1]) * y[3] * y[4];
}

double func_b1(double *y) {
    return y[1] - phi_rad_i;
}
double func_b2(double *y) {
    return y[2] - lambda_rad_i;
}
double func_b3(double *y) {
    return y[1] - phi_rad_f;
}
double func_b4(double *y) {
    return y[2] - lambda_rad_f;
}

double func_dg11(double *y) {
    return 0.;
}
double func_dg12(double *y) {
    return 0.;
}
double func_dg13(double *y) {
    return 1.;
}
double func_dg14(double *y) {
    return 0.;
}
double func_dg21(double *y) {
    return 0.;
}
double func_dg22(double *y) {
    return 0.;
}
double func_dg23(double *y) {
    return 0.;
}
double func_dg24(double *y) {
    return 1.;
}
double func_dg31(double *y) {
    return y[4] * y[4] *
        (sin(y[1]) * sin(y[1]) - cos(y[1]) * cos(y[1]));
}
double func_dg32(double *y) {
    return 0.;
}
double func_dg33(double *y) {
    return 0.;
}
double func_dg34(double *y) {
    return -2. * cos(y[1]) * sin(y[1]) * y[4];
}
double func_dg41(double *y) {
    return 2. * y[3] * y[4] / (cos(y[1]) * cos(y[1]));
}
double func_dg42(double *y) {
    return 0.;
}
double func_dg43(double *y) {
    return 2. * tan(y[1]) * y[4];
}
double func_dg44(double *y) {
    return 2. * tan(y[1]) * y[3];
}

double func_db11(double *y) {return 1.;}
double func_db12(double *y) {return 0.;}
double func_db13(double *y) {return 0.;}
double func_db14(double *y) {return 0.;}
double func_db21(double *y) {return 0.;}
double func_db22(double *y) {return 1.;}
double func_db23(double *y) {return 0.;}
double func_db24(double *y) {return 0.;}
double func_db31(double *y) {return 1.;}
double func_db32(double *y) {return 0.;}
double func_db33(double *y) {return 0.;}
double func_db34(double *y) {return 0.;}
double func_db41(double *y) {return 0.;}
double func_db42(double *y) {return 1.;}
double func_db43(double *y) {return 0.;}
double func_db44(double *y) {return 0.;}

void GeodesicSphere::init(double phi_deg_i, double lambda_deg_i,
                          double phi_deg_f, double lambda_deg_f,
                          double lambda_deg_b,
                          double xi_min_in,
                          double xi_max_in,
                          double delta_xi_in) {
    free_odesolve();

    phi_ini_ = phi_deg_i * fac_deg_to_rad_;
    phi_fin_ = phi_deg_f * fac_deg_to_rad_;

    lambda_ini_ = lambda_deg_i * fac_deg_to_rad_;
    lambda_fin_ = lambda_deg_f * fac_deg_to_rad_;

    lambda_base_ = lambda_deg_b * fac_deg_to_rad_;

    find_lambda_within();

    phi_rad_i = phi_ini_;
    phi_rad_f = phi_fin_;

    lambda_rad_i = lambda_ini_;
    lambda_rad_f = lambda_fin_;

    ptr_odesol_ = new ODESolve::LibRX();

    ptr_odesol_->ptr_log_ = stderr;

    ptr_odesol_->n_func_y_ = 4;
    ptr_odesol_->n_boundary_i_ = 2;

    xi_min_ = xi_min_in;
    xi_max_ = xi_max_in;
    delta_xi_ = delta_xi_in;

    int n_bin_xi_in =
        (int)ceil(fabs(xi_max_ - xi_min_) / delta_xi_);
    delta_xi_ = fabs(xi_max_ - xi_min_) / 
               static_cast<double>(n_bin_xi_in);
    ptr_odesol_->n_bin_x_ = n_bin_xi_in;

    ptr_odesol_->alloc_func();
    for (int ixi = 0; ixi <= ptr_odesol_->n_bin_x_; ixi++) {
        double xi_now = xi_min_ + delta_xi_ * static_cast<double>(ixi);

        ptr_odesol_->tab_func_y_[0][ixi] = xi_now;

        ptr_odesol_->tab_func_y_[1][ixi] =
            (phi_ini_ * (xi_max_ - xi_now) +
             phi_fin_ * (xi_now - xi_min_)) / (xi_max_ - xi_min_);
        ptr_odesol_->tab_func_y_[2][ixi] =
            (lambda_ini_ * (xi_max_ - xi_now) +
             lambda_fin_ * (xi_now - xi_min_)) / (xi_max_ - xi_min_);
        ptr_odesol_->tab_func_y_[3][ixi] =
            (phi_fin_ - phi_ini_) / (xi_max_ - xi_min_);
        ptr_odesol_->tab_func_y_[4][ixi] =
            (lambda_fin_ - lambda_ini_) / (xi_max_ - xi_min_);
    }

    ptr_odesol_->conv_ = 1.0e-5;
    ptr_odesol_->slowc_ = 0.05;

    ptr_odesol_->scalv_[1] = M_PI;
    ptr_odesol_->scalv_[2] = M_PI;
    ptr_odesol_->scalv_[3] = M_PI;
    ptr_odesol_->scalv_[4] = M_PI;

    ptr_odesol_->indexv_[1] = 1;
    ptr_odesol_->indexv_[2] = 2;
    ptr_odesol_->indexv_[3] = 3;
    ptr_odesol_->indexv_[4] = 4;

    ptr_odesol_->ptr_func_g_[1] = &func_g1;
    ptr_odesol_->ptr_func_g_[2] = &func_g2;
    ptr_odesol_->ptr_func_g_[3] = &func_g3;
    ptr_odesol_->ptr_func_g_[4] = &func_g4;

    ptr_odesol_->ptr_func_b_[1] = &func_b1;
    ptr_odesol_->ptr_func_b_[2] = &func_b2;
    ptr_odesol_->ptr_func_b_[3] = &func_b3;
    ptr_odesol_->ptr_func_b_[4] = &func_b4;

    ptr_odesol_->ptr_func_dgdy_[1][1] = &func_dg11;
    ptr_odesol_->ptr_func_dgdy_[1][2] = &func_dg12;
    ptr_odesol_->ptr_func_dgdy_[1][3] = &func_dg13;
    ptr_odesol_->ptr_func_dgdy_[1][4] = &func_dg14;
    ptr_odesol_->ptr_func_dgdy_[2][1] = &func_dg21;
    ptr_odesol_->ptr_func_dgdy_[2][2] = &func_dg22;
    ptr_odesol_->ptr_func_dgdy_[2][3] = &func_dg23;
    ptr_odesol_->ptr_func_dgdy_[2][4] = &func_dg24;
    ptr_odesol_->ptr_func_dgdy_[3][1] = &func_dg31;
    ptr_odesol_->ptr_func_dgdy_[3][2] = &func_dg32;
    ptr_odesol_->ptr_func_dgdy_[3][3] = &func_dg33;
    ptr_odesol_->ptr_func_dgdy_[3][4] = &func_dg34;
    ptr_odesol_->ptr_func_dgdy_[4][1] = &func_dg41;
    ptr_odesol_->ptr_func_dgdy_[4][2] = &func_dg42;
    ptr_odesol_->ptr_func_dgdy_[4][3] = &func_dg43;
    ptr_odesol_->ptr_func_dgdy_[4][4] = &func_dg44;

    ptr_odesol_->ptr_func_dbdy_[1][1] = &func_db11;
    ptr_odesol_->ptr_func_dbdy_[1][2] = &func_db12;
    ptr_odesol_->ptr_func_dbdy_[1][3] = &func_db13;
    ptr_odesol_->ptr_func_dbdy_[1][4] = &func_db14;
    ptr_odesol_->ptr_func_dbdy_[2][1] = &func_db21;
    ptr_odesol_->ptr_func_dbdy_[2][2] = &func_db22;
    ptr_odesol_->ptr_func_dbdy_[2][3] = &func_db23;
    ptr_odesol_->ptr_func_dbdy_[2][4] = &func_db24;
    ptr_odesol_->ptr_func_dbdy_[3][1] = &func_db31;
    ptr_odesol_->ptr_func_dbdy_[3][2] = &func_db32;
    ptr_odesol_->ptr_func_dbdy_[3][3] = &func_db33;
    ptr_odesol_->ptr_func_dbdy_[3][4] = &func_db34;
    ptr_odesol_->ptr_func_dbdy_[4][1] = &func_db41;
    ptr_odesol_->ptr_func_dbdy_[4][2] = &func_db42;
    ptr_odesol_->ptr_func_dbdy_[4][3] = &func_db43;
    ptr_odesol_->ptr_func_dbdy_[4][4] = &func_db44;

    ptr_odesol_->init();

    int n_size_tab_xi = ptr_odesol_->n_bin_x_ + 1;
    tab_xi_.resize(n_size_tab_xi);
    tab_phi_.resize(n_size_tab_xi);
    tab_lambda_.resize(n_size_tab_xi);
    tab_distance_.resize(n_size_tab_xi);

    initialized_ = true;

    return;
}

void GeodesicSphere::find_geodesic() {
    if (!initialized_) {
        return;
    }

    bool converging = false;
    while (!converging) {
        converging = ptr_odesol_->next();
    }

    double **xvec = new double *[2];
    xvec[0] = new double[3];
    xvec[1] = new double[3];

    for (int ixi = 0; ixi <= ptr_odesol_->n_bin_x_; ixi++) {
        tab_xi_[ixi] = ptr_odesol_->tab_func_y_[0][ixi];

        tab_phi_[ixi] = ptr_odesol_->tab_func_y_[1][ixi];
        tab_lambda_[ixi] = ptr_odesol_->tab_func_y_[2][ixi];

        if (ixi == 0) {
            tab_distance_[ixi] = 0.;
        } else {
            for (int jxi = 0; jxi < 2; jxi++) {
                xvec[jxi][0] = sin(tab_phi_[ixi - jxi]);
                xvec[jxi][1] = cos(tab_phi_[ixi - jxi]) *
                               cos(tab_lambda_[ixi - jxi]);
                xvec[jxi][2] = cos(tab_phi_[ixi - jxi]) *
                               sin(tab_lambda_[ixi - jxi]);
            }

            double cos_dist =
                xvec[0][0] * xvec[1][0] +
                xvec[0][1] * xvec[1][1] +
                xvec[0][2] * xvec[1][2];

            tab_distance_[ixi] =
                tab_distance_[ixi - 1] + acos(cos_dist);
        }
    }

    delete [] xvec[0];
    delete [] xvec[1];
    delete [] xvec;

    return;
}

void GeodesicSphere::find_lambda_within() {
    bool found_lambda_i =
        fabs(lambda_ini_ - lambda_base_) < M_PI;
    while (!found_lambda_i) {
        if (lambda_ini_ < lambda_base_) {
            lambda_ini_ += 2. * M_PI;
        } else {
            lambda_ini_ -= 2. * M_PI;
        }

        found_lambda_i =
            fabs(lambda_ini_ - lambda_base_) < M_PI;
    }

    bool found_lambda_f =
        fabs(lambda_fin_ - lambda_ini_) < M_PI;
    while (!found_lambda_f) {
        if (lambda_fin_ < lambda_ini_) {
            lambda_fin_ += 2. * M_PI;
        } else {
            lambda_fin_ -= 2. * M_PI;
        }

        found_lambda_f =
            fabs(lambda_fin_ - lambda_ini_) < M_PI;
    }

    return;
}
