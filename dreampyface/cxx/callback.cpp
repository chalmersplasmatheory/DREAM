/**
 * Handling of callbacks from libDREAM.
 */

#include "pyface/callback.hpp"
#include "DREAM/Simulation.hpp"

std::vector<PyObject*>
    callback_timestep,
    callback_iteration,
    callback_timestep_term;


/**
 * Put the 'DREAM::Simulation' object into a capsule and
 * construct a new 'dreampyface.Simulation' object.
 */
PyObject *capsule_to_simulation(DREAM::Simulation *sim) {
    PyObject *cap = PyCapsule_New(sim, "sim", NULL);

    //PyObject *name = PyUnicode_FromString("Simulation");
	const char *name = "Simulation";
    PyObject *sys_mod_dict = PyImport_GetModuleDict();

    PyObject *sim_mod = PyMapping_GetItemString(sys_mod_dict, "dreampyface.Simulation");
    if (sim_mod == nullptr) {
        PyErr_PrintEx(1);
        throw DREAM::DREAMException(
            "Python termination function error."
        );
    }

    //PyObject *pysim = PyObject_CallMethodOneArg(sim_mod, name, cap);
	PyObject *pysim = PyObject_CallMethod(sim_mod, name, "O", cap);
    if (pysim == nullptr) {
        PyErr_PrintEx(1);
        throw DREAM::DREAMException(
            "Python termination function error."
        );
    }

    Py_DECREF(cap);

    return pysim;
}

/**
 * Helper function for registering the DREAMpy callback functions
 * on the given DREAM::Simulation object.
 */
void register_callback_functions(DREAM::Simulation *sim) {
    DREAM::EquationSystem *eqsys = sim->GetEquationSystem();

    if (callback_timestep.size() > 0)
        eqsys->RegisterCallback_TimestepFinished(&dreampy_callback_timestep);

    if (callback_iteration.size() > 0)
        eqsys->RegisterCallback_IterationFinished(&dreampy_callback_iteration);
}

/**
 * Central callback function for when DREAM has advanced the
 * equation system by one time step.
 */
void dreampy_callback_timestep(DREAM::Simulation *sim) {
    PyObject *pysim = capsule_to_simulation(sim);
	PyObject *args = PyTuple_Pack(1, pysim);

    for (auto f : callback_timestep)
		PyObject_Call(f, args, NULL);

	Py_DECREF(args);
    Py_DECREF(pysim);
}

/**
 * Central callback function for when the DREAM non-linear solver
 * has completed another iteration.
 */
void dreampy_callback_iteration(DREAM::Simulation *sim) {
    PyObject *pysim = capsule_to_simulation(sim);
	PyObject *args = PyTuple_Pack(1, pysim);

    for (auto f : callback_iteration)
        PyObject_Call(f, args, NULL);

	Py_DECREF(args);
    Py_DECREF(pysim);
}

/**
 * Central callback function for determining when the time
 * stepping should be terminated.
 */
bool dreampy_callback_return_bool(void *func, DREAM::Simulation *sim) {
    PyObject *pysim = capsule_to_simulation(sim);
	PyObject *args = PyTuple_Pack(1, pysim);
    bool v = true;

    PyObject *f = (PyObject*)func;

    PyObject *ret = PyObject_Call(f, args, NULL);

    if (ret == nullptr) {
        PyErr_PrintEx(1);
        throw DREAM::DREAMException(
            "Python termination function error."
        );
    }

    v = PyObject_IsTrue(ret);
    Py_DECREF(ret);

	Py_DECREF(args);
    Py_DECREF(pysim);

    return v;
}

