#ifndef _WRAPGEOIATA_H_
#define _WRAPGEOIATA_H_

#include<string>
#include<Python.h>

class WrapGeoIATA {
  private :

    PyObject *ptr_py_instance_ = NULL;

    std::string name_city_;
    std::string name_country_;
    double lat_deg_;
    double lon_deg_;

    bool initialized_;

  public :

    std::string path_python_;
    std::string path_module_;

    WrapGeoIATA() {
        path_python_ = "";
        path_module_ = "";

        name_city_ = "";
        name_country_ = "";
        lat_deg_ = 0.;
        lon_deg_ = 0.;

        initialized_ = false;

        return;
    }

    ~WrapGeoIATA() {
        if (!initialized_) {
            return;
        }

        Py_Finalize();

        return;
    }

    void init();

    void set_location(std::string &code_iata);

    void verbose(FILE *ptr_log = stdout);

    void get_city(std::string *ptr_name) {
        *ptr_name = name_city_;

        return;
    }

    void get_country(std::string *ptr_name) {
        *ptr_name = name_country_;

        return;
    }

    double get_latitude() {return lat_deg_;}
    double get_longitude() {return lon_deg_;}
};

#endif
