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

    /* tabulated t / tau,
     * where tau = sqrt(pow(a, 3) / G / M)
     *   a = semi-major axis
     *   M = mass of the Sun (or center of the gravitational force)
     *   G = Gravitational constant
     *
     * Note that the period of the orbit is T = 2 pi tau. */
    std::vector<double> tab_t_;
    /* tabulated x
     * divided by the semi-major axis */
    std::vector<double> tab_x_;
    /* tabulated y
     * divided by the semi-major axis */
    std::vector<double> tab_y_;

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

        tab_t_.clear();
        tab_x_.clear();
        tab_y_.clear();

        ptr_odesol_->free_func();
        delete ptr_odesol_;

        initialized_ = false;

        return;
    }

    /* initializing function,
     * which must be called first
     *
     * The following arguments must be provided.
     *   r_i_in : initial distance from the center
     *            divided by the semi-major axis
     *   vr_i_in : initial radial velocity
     *   phi_i_in : initial azimuthal angle on the xy-plane
     *   l_ang_in : angular momentum per unit mass
     *              divided by sqrt(G M a)
     *   t_min_in : initial time
     *              divided by tau = sqrt(pow(a, 3) / G / M)
     *       Default value is 0.
     *   t_max_in : final time
     *              divided by tau = sqrt(pow(a, 3) / G / M)
     *       Default value is 2 pi
     *   delta_t_in : increment in time for each Runge-Kutta step
     *       Default value is 0.01 */
    void init(double r_i_in,
              double vr_i_in,
              double phi_i_in,
              double l_ang_in,
              double t_min_in = 0.,
              double t_max_in = 2. * M_PI,
              double delta_t_in = 0.01);

    /* solve differential equations
     * to obtain the orbit
     * and populate tabulated functions */
    void find_orbit();

    /* write out an orbit on text file
     *
     * The output file contains
     * x / a and y / a as functions of time / tau. */
    void export_orbit(char *filename) {
        if (!initialized_) {
            return;
        }

        FILE *ptr_fout;
        ptr_fout = fopen(filename, "w");
        if (ptr_fout == NULL) {
            return;
        }

        fprintf(ptr_fout, "# output from OrbitKepler\n");
        fprintf(ptr_fout, "# t/tau    ");
        fprintf(ptr_fout, "x/a    ");
        fprintf(ptr_fout, "y/a\n");

        for (int it = 0; it <= n_bin_t_; it++) {
            fprintf(ptr_fout, "    %e    %e    %e\n",
                    tab_t_[it], tab_x_[it], tab_y_[it]);
        }

        fclose(ptr_fout);

        return;
    }
};

#endif
