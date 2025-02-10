#ifndef _GEODESICSPHERE_H_
#define _GEODESICSPHERE_H_

#include<math.h>
#include"ODESolveRX.h"

class GeodesicSphere {
  private :

    double fac_deg_to_rad_;

    double phi_ini_;
    double phi_fin_;

    double lambda_ini_;
    double lambda_fin_;

    double lambda_base_;

    double xmin_;
    double xmax_;
    double delta_x_;

    ODESolve::LibRX *ptr_odesol_;

    bool initialized_;

  public :

    GeodesicSphere() {
        fac_deg_to_rad_ = M_PI / 180.;

        initialized_ = false;

        return;
    }

    ~GeodesicSphere() {
        free_odesolve();

        return;
    }

    void free_odesolve() {
        if (!initialized_) {
            return;
        }

        ptr_odesol_->free_func();
        delete ptr_odesol_;

        initialized_ = false;

        return;
    }

    void init(double phi_deg_i, double lambda_deg_i,
              double phi_deg_f, double lambda_deg_f,
              double lambda_deg_b = 0.,
              double xmin_in = 0.,
              double xmax_in = 1.,
              double delta_x_in = 0.0005);

    void find_geodesic();

    void export_geodesic(char *filename) {
        if (!initialized_) {
            return;
        }

        FILE *ptr_fout;
        ptr_fout = fopen(filename, "w");
        if (ptr_fout == NULL) {
            return;
        }

        ptr_odesol_->export_file(ptr_fout);

        fclose(ptr_fout);

        return;
    }

    void find_lambda_within();
};

#endif