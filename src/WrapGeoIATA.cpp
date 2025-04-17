#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"PythonHook.h"
#include"WrapGeoIATA.h"

void WrapGeoIATA::init() {
    std::string name_module = "GeoIATA";
    PyObject *ptr_py_dict =
        PythonHook::get_ptr_dict(name_module);

    std::string name_class = "GeoLocation";
    PyObject *ptr_py_class =
        PythonHook::get_ptr_class(ptr_py_dict, name_class);

    ptr_py_GeoIATA_ =
        PythonHook::get_ptr_instance(ptr_py_class, NULL);
    if (ptr_py_GeoIATA_ == NULL) {
        return;
    }

    initialized_ = true;

    return;
}

void WrapGeoIATA::set_location(std::string &code_iata) {
    if (!initialized_ || ptr_py_GeoIATA_ == NULL) {
        return;
    }

    PyObject *p_value_input;
    PyObject *p_value_found;
    PyObject *p_value_city;
    PyObject *p_value_country;
    PyObject *p_value_lat;
    PyObject *p_value_lon;

    p_value_input =
        PyObject_CallMethod(ptr_py_GeoIATA_,
                            "import_airport",
                            "s", code_iata.c_str());

    p_value_found =
        PyObject_CallMethod(ptr_py_GeoIATA_,
                            "get_found", NULL);
    is_found_ = PyObject_IsTrue(p_value_found) != 0;

    p_value_city =
        PyObject_CallMethod(ptr_py_GeoIATA_,
                            "get_city", NULL);
    name_city_ = PyUnicode_AsUTF8(p_value_city);

    p_value_country =
        PyObject_CallMethod(ptr_py_GeoIATA_,
                            "get_country", NULL);
    name_country_ = PyUnicode_AsUTF8(p_value_country);

    p_value_lat = PyObject_CallMethod(ptr_py_GeoIATA_,
                                      "get_latitude", NULL);
    lat_deg_ = PyFloat_AS_DOUBLE(p_value_lat);

    p_value_lon = PyObject_CallMethod(ptr_py_GeoIATA_,
                                      "get_longitude", NULL);
    lon_deg_ = PyFloat_AS_DOUBLE(p_value_lon);

    return;
}

void WrapGeoIATA::verbose(FILE *ptr_log) {
    if (!initialized_ || ptr_log == NULL) {
        return;
    }

    fprintf(ptr_log, "  %s, %s\n",
        name_city_.c_str(), name_country_.c_str());
    fprintf(ptr_log, "    lat (deg) = %e\n", lat_deg_);
    fprintf(ptr_log, "    lon (deg) = %e\n", lon_deg_);

    return;
}
