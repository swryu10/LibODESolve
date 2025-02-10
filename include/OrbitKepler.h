#ifndef _ORBITKEPLER_H_
#define _ORBITKEPLER_H_

#include<math.h>
#include"ODESolveRK.h"

class OrbitKepler {
  private :

    double r_ini_;
    double vr_ini_;
    double phi_ini_;

    double l_angular_;

    int n_bin_t_;
    double t_min_;
    double t_max_;
    double delta_t_;

    ODESolve::LibRK *ptr_odesol_;

    bool initialized_;

  public :

    OrbitKepler() {
        initialized_ = false;

        return;
    }

    ~OrbitKepler() {
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

    void init(double r_i_in,
              double vr_i_in,
              double phi_i_in,
              double l_ang_in,
              double t_min_in = 0.,
              double t_max_in = 2. * M_PI,
              double delta_t_in = 0.01);

    void find_orbit();

    void export_orbit(char *filename) {
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
};

#endif
