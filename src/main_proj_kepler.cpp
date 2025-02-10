#include"OrbitKepler.h"

int main(int argc, char *argv[]) {
    OrbitKepler orbit_kepler;

    double ecc = 0.6;

    double r_i = 1. - ecc;
    double vr_i = 0.;
    double phi_i = 0.;

    double l_ang = sqrt(1. - ecc * ecc);

    double t_min = 0.;
    double t_max = 2. * M_PI;
    double delta_t = 0.01;

    orbit_kepler.init(r_i, vr_i, phi_i, l_ang,
                      t_min, t_max, delta_t);

    orbit_kepler.find_orbit();

    char fname_out[200];
    strcpy(fname_out, "orbit_kepler.txt");

    orbit_kepler.export_orbit(fname_out);

    return 0;
}
