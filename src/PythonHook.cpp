#include"PythonHook.h"

namespace PythonHook {

std::string path_python_ = "./";
std::string path_module_ = "./";

void func_py_ini() {
    setenv("PYTHONPATH", path_python_.c_str(), 1);

    char state_sys_append[100];
    strcpy(state_sys_append, "sys.path.append(\'");
    strcat(state_sys_append, path_module_.c_str());
    strcat(state_sys_append, "\')\n");

    Py_Initialize();
    PyRun_SimpleString("import sys\n");
    PyRun_SimpleString(state_sys_append);

    return;
}

void func_py_fin() {
    Py_Finalize();

    return;
}

PyObject *get_ptr_dict(std::string name_module) {
    PyObject *ptr_py_name =
        PyUnicode_DecodeFSDefault(name_module.c_str());
    if (ptr_py_name == NULL) {
        return NULL;
    }

    PyObject *ptr_py_module =
        PyImport_Import(ptr_py_name);
    if (ptr_py_module == NULL) {
        return NULL;
    }

    Py_DECREF(ptr_py_name);

    PyObject *ptr_py_dict =
        PyModule_GetDict(ptr_py_module);

    Py_DECREF(ptr_py_module);

    return ptr_py_dict;
}

PyObject *get_ptr_class(PyObject *ptr_py_dict,
                        std::string name_class) {
    PyObject *ptr_py_class =
        PyDict_GetItemString(ptr_py_dict, name_class.c_str());

    if (PyCallable_Check(ptr_py_class)) {
        return ptr_py_class;
    } else {
        return NULL;
    }
}

PyObject *get_ptr_instance(PyObject *ptr_py_class,
                           PyObject *ptr_py_args) {
    if (ptr_py_class == NULL) {
        return NULL;
    }

    PyObject *ptr_py_instance =
        PyObject_CallObject(ptr_py_class, ptr_py_args);

    return ptr_py_instance;
}

} // end namespace PythonHook
