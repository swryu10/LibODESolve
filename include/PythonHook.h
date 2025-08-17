#ifndef _PYTHONHOOK_H_
#define _PYTHONHOOK_H_

#include<string>
#include<Python.h>

class PythonHook {
  private :

    static std::string path_python_;
    static std::string path_module_;

  public :

    static void func_py_ini();
    static void func_py_fin();

    static void set_path_python(std::string path_in);
    static void set_path_module(std::string path_in);

    static PyObject *get_ptr_dict(std::string name_module);
    static PyObject *get_ptr_class(PyObject *ptr_py_dict,
                                   std::string name_class);
    static PyObject *get_ptr_instance(PyObject *ptr_py_class,
                                      PyObject *ptr_py_args = NULL);

};

#endif
