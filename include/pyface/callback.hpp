#ifndef _DREAM_PYFACE_CALLBACK_HPP
#define _DREAM_PYFACE_CALLBACK_HPP

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include <vector>
#include "DREAM/Simulation.hpp"

extern std::vector<PyObject*>
    callback_timestep,
    callback_iteration;

void register_callback_functions(DREAM::Simulation*);
PyObject *capsule_to_simulation(DREAM::Simulation*);
void dreampy_callback_timestep(DREAM::Simulation*);
void dreampy_callback_iteration(DREAM::Simulation*);
bool dreampy_callback_return_bool(void*, DREAM::Simulation*);

#endif/*_DREAM_PYFACE_CALLBACK_HPP*/
