#ifndef _GEODESICSPHERE_H_
#define _GEODESICSPHERE_H_

#include<math.h>
#include<vector>
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

    /* tabulated x
     * (parametric variable) */
    std::vector<double> tab_x_;
    /* tabulated latitude
     * in radian */
    std::vector<double> tab_phi_;
    /* tabulated longitude
     * in radian */
    std::vector<double> tab_lambda_;
    /* tabulated distance from the origin
     * divided by radius of sphere */
    std::vector<double> tab_distance_;

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

        tab_x_.clear();
        tab_phi_.clear();
        tab_lambda_.clear();
        tab_distance_.clear();

        ptr_odesol_->free_func();
        delete ptr_odesol_;

        initialized_ = false;

        return;
    }

    /* initializing function,
     * which must be called first
     *
     * The following arguments must be provided.
     *   phi_deg_i : latitude of the origin in deg
     *   lambda_deg_i : longitude of the origin in deg
     *   phi_deg_f : latitude of the destination in deg
     *   lambda_deg_f : longitude of the destination in deg
     *   lambda_deg_b : base longitude
     *       integer multiple of 2 pi is added to or subtracted from
     *       lambda_deg_i, such that
     *       |lambda_deg_i - lambda_deg_b| < pi.
     *       Default value is 0.
     *   xmin_in : minimum value of the parametric variable x
     *       Default value is 0.
     *   xmax_in : maximum value of the parametric variable x
     *       Default value is 1.
     *   delta_x_in : increment in x for each mesh point
     *       Default value is 0.0005 */
    void init(double phi_deg_i, double lambda_deg_i,
              double phi_deg_f, double lambda_deg_f,
              double lambda_deg_b = 0.,
              double xmin_in = 0.,
              double xmax_in = 1.,
              double delta_x_in = 0.0005);

    /* solve differential equations
     * to obtain the geodesic
     * and populate tabulated functions */
    void find_geodesic();

    /* write out a geodesic on text file
     *
     * The output file contains
     * latitude (in deg), longitude (in deg)
     * and distance / radius
     * as functions of the parametric variable x. */
    void export_geodesic(char *filename) {
        if (!initialized_) {
            return;
        }

        FILE *ptr_fout;
        ptr_fout = fopen(filename, "w");
        if (ptr_fout == NULL) {
            return;
        }

        fprintf(ptr_fout, "# output from GeodesicSphere\n");
        fprintf(ptr_fout, "# x    ");
        fprintf(ptr_fout, "latitude(deg)    ");
        fprintf(ptr_fout, "longitude(deg)    ");
        fprintf(ptr_fout, "distance/radius\n");

        for (int ix = 0; ix <= ptr_odesol_->n_bin_x_; ix++) {
            fprintf(ptr_fout, "    %e    %e    %e    %e\n",
                    tab_x_[ix],
                    tab_phi_[ix] / fac_deg_to_rad_,
                    tab_lambda_[ix] / fac_deg_to_rad_,
                    tab_distance_[ix]);
        }

        fclose(ptr_fout);

        return;
    }

    void find_lambda_within();
};

#endif