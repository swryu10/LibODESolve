#ifndef _WRAPGEOIATA_H_
#define _WRAPGEOIATA_H_

#include<string>
#include<Python.h>

class WrapGeoIATA {
  private :

    PyObject *ptr_py_dict_;
    PyObject *ptr_py_class_;
    PyObject *ptr_py_GeoIATA_;

    bool is_found_;
    std::string name_city_;
    std::string name_country_;
    double lat_deg_;
    double lon_deg_;

    bool initialized_;

  public :

    WrapGeoIATA() {
        ptr_py_dict_ = NULL;
        ptr_py_class_ = NULL;
        ptr_py_GeoIATA_ = NULL;

        is_found_ = false;
        name_city_ = "";
        name_country_ = "";
        lat_deg_ = 0.;
        lon_deg_ = 0.;

        initialized_ = false;

        return;
    }

    ~WrapGeoIATA() {}

    void free_ptr_py() {
        if (!initialized_) {
            return;
        }

        Py_XDECREF(ptr_py_dict_);
        Py_XDECREF(ptr_py_class_);
        Py_XDECREF(ptr_py_GeoIATA_);

        return;
    }

    void init();

    void set_location(std::string &code_iata);

    void verbose(FILE *ptr_log = stdout);

    bool get_found() {return is_found_;}

    std::string get_city() {
        return name_city_;
    }

    std::string get_country() {
        return name_country_;
    }

    double get_latitude() {return lat_deg_;}
    double get_longitude() {return lon_deg_;}
};

#endif
