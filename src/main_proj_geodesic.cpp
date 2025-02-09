#include"GeodesicSphere.h"

// latitude (degree)
double phi_deg_ICN = 37.4602;
double phi_deg_JFK = 40.6413;
double phi_deg_LAX = 33.9416;
double phi_deg_HNL = 21.3069;
double phi_deg_SIN = 1.3521;
double phi_deg_LHR = 51.4700;
double phi_deg_FRA = 50.0379;
double phi_deg_SYD = -33.8688;
double phi_deg_GRU = -23.4306;
double phi_deg_JNB = -26.2041;
double phi_deg_DXB = 25.2048;

// longitude (degree)
double lambda_deg_ICN = 126.4407;
double lambda_deg_JFK = (360. - 73.7781);
double lambda_deg_LAX = (360. - 118.4085);
double lambda_deg_HNL = (360. - 157.8583);
double lambda_deg_SIN = 103.8198;
double lambda_deg_LHR = -0.4543;
double lambda_deg_FRA = 8.5622;
double lambda_deg_SYD = 151.2093;
double lambda_deg_GRU = - 46.4730;
double lambda_deg_JNB = 28.0473;
double lambda_deg_DXB = 55.2708;

int main(int argc, char *argv[]) {
    GeodesicSphere geosphere;

    char name_org[] = "ICN";
    char name_dst[] = "LAX";

    double phi_deg_i = phi_deg_ICN;
    double phi_deg_f = phi_deg_LAX;

    double lambda_deg_i = 0.;
    double lambda_deg_f = lambda_deg_LAX - lambda_deg_ICN;

    double lambda_deg_b = 0.;

    double xmin = 0.;
    double xmax = 1.;
    double deltax = 0.0001;

    geosphere.init(phi_deg_i, lambda_deg_i,
                   phi_deg_f, lambda_deg_f,
                   lambda_deg_b,
                   xmin, xmax, deltax);
    geosphere.find_geodesic();

    char fname_out[200];
    strcpy(fname_out, "geodesic_");
    strcat(fname_out, name_org);
    strcat(fname_out, "_");
    strcat(fname_out, name_dst);
    strcat(fname_out, ".txt");

    geosphere.export_geodesic(fname_out);

    return 0;
}
