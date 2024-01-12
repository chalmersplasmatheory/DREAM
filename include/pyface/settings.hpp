#ifndef _DREAM_PYFACE_SETTINGS_HPP
#define _DREAM_PYFACE_SETTINGS_HPP

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include "DREAM/Settings/Settings.hpp"

DREAM::Settings *dreampy_loadsettings(PyObject*);
void dreampy_load_dict(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_address(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_bool(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_int(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_real(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_int_array(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_real_array(DREAM::Settings*, const std::string&, PyObject*);
void dreampy_load_string(DREAM::Settings*, const std::string&, PyObject*);

#endif/*_DREAM_PYFACE_SETTINGS_HPP*/
