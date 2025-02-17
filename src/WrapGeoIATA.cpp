#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"WrapGeoIATA.h"

void WrapGeoIATA::init() {
    if (initialized_) {
        return;
    }

    setenv("PYTHONPATH", path_python_.c_str(), 1);

    char state_sys_append[100];
    strcpy(state_sys_append, "sys.path.append(\'");
    strcat(state_sys_append, path_module_.c_str());
    strcat(state_sys_append, "\')\n");

    Py_Initialize();
    PyRun_SimpleString("import sys\n");
    PyRun_SimpleString(state_sys_append);

    PyObject *ptr_py_name = NULL;
    PyObject *ptr_py_module = NULL;
    PyObject *ptr_py_dict = NULL;
    PyObject *ptr_py_class = NULL;

    ptr_py_name = PyUnicode_DecodeFSDefault("GeoIATA");
    if (ptr_py_name == NULL) {
        return;
    }

    ptr_py_module = PyImport_Import(ptr_py_name);
    if (ptr_py_module == NULL) {
        return;
    }

    Py_DECREF(ptr_py_name);

    ptr_py_dict = PyModule_GetDict(ptr_py_module);
    if (ptr_py_dict == NULL) {
        return;
    }

    Py_DECREF(ptr_py_module);

    ptr_py_class =
        PyDict_GetItemString(ptr_py_dict, "GeoLocation");
    if (PyCallable_Check(ptr_py_class)) {
        ptr_py_instance_ =
            PyObject_CallObject(ptr_py_class, NULL);
    }

    initialized_ = true;

    return;
}

void WrapGeoIATA::set_location(char *code_iata) {
    if (!initialized_ || ptr_py_instance_ == NULL) {
        return;
    }

    PyObject *p_value_input;
    PyObject *p_value_city;
    PyObject *p_value_country;
    PyObject *p_value_lat;
    PyObject *p_value_lon;

    p_value_input =
        PyObject_CallMethod(ptr_py_instance_,
                            "import_airport",
                            "s", code_iata);

    p_value_city =
        PyObject_CallMethod(ptr_py_instance_,
                            "get_city", NULL);
    name_city_ = PyUnicode_AsUTF8(p_value_city);

    p_value_country =
        PyObject_CallMethod(ptr_py_instance_,
                            "get_country", NULL);
    name_country_ = PyUnicode_AsUTF8(p_value_country);

    p_value_lat = PyObject_CallMethod(ptr_py_instance_,
                                      "get_latitude", NULL);
    lat_deg_ = PyFloat_AS_DOUBLE(p_value_lat);

    p_value_lon = PyObject_CallMethod(ptr_py_instance_,
                                      "get_longitude", NULL);
    lon_deg_ = PyFloat_AS_DOUBLE(p_value_lon);

    return;
}
