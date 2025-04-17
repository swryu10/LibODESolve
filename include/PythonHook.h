#ifndef _PYTHONHOOK_H_
#define _PYTHONHOOK_H_

#include<string>
#include<Python.h>

namespace PythonHook {

extern std::string path_python_;
extern std::string path_module_;

void func_py_ini();
void func_py_fin();

PyObject *get_ptr_dict(std::string name_module);
PyObject *get_ptr_class(PyObject *ptr_py_dict,
                        std::string name_class);
PyObject *get_ptr_instance(PyObject *ptr_py_class,
                           PyObject *ptr_py_args = NULL);

} // end namespace PythonHook

#endif
