#include"GeodesicSphere.h"
#include"WrapGeoIATA.h"

int main(int argc, char *argv[]) {
    GeodesicSphere geosphere;
    WrapGeoIATA geolocation;

    std::string name_org;
    if (argc > 1) {
        name_org = argv[1];
    } else {
        char name_in[100];
        fprintf(stdout, "IATA code for the origin : ");
        scanf("%s", name_in);

        name_org = name_in;
    }

    std::string name_dst;
    if (argc > 2) {
        name_dst = argv[2];
    } else {
        char name_in[100];
        fprintf(stdout, "IATA code for the destination : ");
        scanf("%s", name_in);

        name_dst = name_in;
    }

    fprintf(stdout, "\n");

    char env_libode_path_python[1000];
    strcpy(env_libode_path_python,
           getenv("LIBODE_PATH_PYTHON"));
    //fprintf(stderr, "  %s\n", env_libode_path_python);
    geolocation.path_python_ = env_libode_path_python;
    geolocation.path_module_ = "./";

    geolocation.init();

    fprintf(stdout, "Origin :\n");
    geolocation.set_location(name_org);
    if (!geolocation.get_found()) {
        fprintf(stderr, "  IATA code %s not found.\n",
                name_org.c_str());
        fprintf(stderr, "  Program aborted.\n");

        return 0;
    }
    geolocation.verbose();
    double phi_deg_i = geolocation.get_latitude();
    double lambda_deg_i = geolocation.get_longitude();

    fprintf(stdout, "Destination :\n");
    geolocation.set_location(name_dst);
    if (!geolocation.get_found()) {
        fprintf(stderr, "  IATA code %s not found.\n",
                name_dst.c_str());

        fprintf(stderr, "  Program aborted.\n");

        return 0;
    }
    geolocation.verbose();
    double phi_deg_f = geolocation.get_latitude();
    double lambda_deg_f = geolocation.get_longitude();

    double lambda_deg_b = lambda_deg_i;

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
    strcat(fname_out, name_org.c_str());
    strcat(fname_out, "_");
    strcat(fname_out, name_dst.c_str());
    strcat(fname_out, ".txt");

    geosphere.export_geodesic(fname_out);

    return 0;
}
