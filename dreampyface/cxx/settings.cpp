/**
 * This module contains routines for converting a Python
 * dictionary into a DREAM::Settings object.
 */

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>

/*#ifndef NPY_NO_DEPRECATED_API
#   define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif*/
#include "pyface/numpy.h"
#include <iostream>
#include "DREAM/DREAMException.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "pyface/settings.hpp"
#include "pyface/dreamtypes.hpp"


using DREAM::Settings;
using namespace std;


/**
 * Creates a new 'DREAM::Settings' object from the given
 * Python dict object.
 *
 * dict: Object representing a Python dictionary which contains
 *       the settings to populate the 'Settings' object with.
 */
DREAM::Settings *dreampy_loadsettings(PyObject *dict) {
    DREAM::Settings *set = DREAM::SimulationGenerator::CreateSettings();

    dreampy_load_dict(set, "", dict);

    return set;
}

/**
 * Function which is used to recursively traverse and parse
 * the Python dictionary containing the settings.
 * 
 * s:    Settings object to populate.
 * path: Path in settings currently being traversed.
 * dict: Python dictionary to parse.
 */
void dreampy_load_dict(DREAM::Settings *s, const string& path, PyObject *dict) {
    PyObject *keys   = PyDict_Keys(dict);
    Py_ssize_t nkeys = PyList_Size(keys);

    for (Py_ssize_t i = 0; i < nkeys; i++) {
        // Get key
        PyObject *key = PyList_GetItem(keys, i);
        const char *keyname = PyUnicode_AsUTF8(key);

        // Construct full name of setting
        string sname;
        if (path.size() == 0)
            sname = keyname;
        else
            sname = path + '/' + keyname;

        //cout << "Loading '" << sname << "'..." << endl;

        // Get dictionary value
        PyObject *val = PyDict_GetItem(dict, key);

        // Check if this item contains data directly, or
        // if it represents another dictionary
        if (PyDict_Check(val)) {
            dreampy_load_dict(s, sname, val);
        } else {
            if (!s->HasSetting(sname))
                continue;

            // Item is a value...
            enum Settings::setting_type tp = s->GetType(sname);

            switch (tp) {
                case Settings::SETTING_TYPE_BOOL:       dreampy_load_bool(s, sname, val); break;
                case Settings::SETTING_TYPE_INT:        dreampy_load_int(s, sname, val); break;
                case Settings::SETTING_TYPE_REAL:       dreampy_load_real(s, sname, val); break;
                case Settings::SETTING_TYPE_STRING:     dreampy_load_string(s, sname, val); break;
                case Settings::SETTING_TYPE_INT_ARRAY:  dreampy_load_int_array(s, sname, val); break;
                case Settings::SETTING_TYPE_REAL_ARRAY: dreampy_load_real_array(s, sname, val); break;
                case Settings::SETTING_TYPE_ADDRESS:    dreampy_load_address(s, sname, val); break;

                default:
                    throw DREAM::DREAMException(
                        "Setting '%s': Unrecognized setting type: %d.",
                        sname.c_str(), tp
                    );
            }
        }
    }
}

/**
 * Load a setting as an address from the given Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_address(Settings *s, const string& name, PyObject *obj) {
    if (PyFunction_Check(obj)) {
        Py_INCREF(obj);
        s->SetSetting(name, (void*)obj);
    } else
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value: %s. Expected function.",
            name.c_str(), obj->ob_type->tp_name
        );
}

/**
 * Load a setting as a bool from the given
 * Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_bool(Settings *s, const string& name, PyObject *obj) {
    if (PyBool_Check(obj)) {
        s->SetSetting(name, reinterpret_cast<bool>(obj == Py_True));
    } else if (PyLong_Check(obj)) {
        long l = PyLong_AsLong(obj);
        s->SetSetting(name, reinterpret_cast<bool>(l != 0));
    } else if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        int ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        if (ndim != 1 || _dims[0] != 1)
            throw DREAM::DREAMException(
                "Setting '%s': Expected value to be a scalar, but array was given.",
                name.c_str()
            );

        int dtype = PyArray_TYPE(ao);
        int_t v;
        if (dtype == NPY_INT) {
            v = reinterpret_cast<int*>(PyArray_DATA(ao))[0];
        } else if (dtype == NPY_LONG) {
            v = reinterpret_cast<long*>(PyArray_DATA(ao))[0];
        } else
            throw DREAM::DREAMException(
                "Setting '%s': Unrecognized data type of specified value: %d. Expected numpy integer array. (1)",
                name.c_str(), dtype
            );

        s->SetSetting(name, v!=0);
    } else
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value: %s. Expected boolean or integer.",
            name.c_str(), obj->ob_type->tp_name
        );
}

/**
 * Load a setting as an integer from the
 * given Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_int(Settings *s, const string& name, PyObject *obj) {
    if (PyLong_Check(obj)) {
        long l = PyLong_AsLong(obj);
        s->SetSetting(name, static_cast<int_t>(l));
    } else if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        int ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        if (ndim != 1 || _dims[0] != 1)
            throw DREAM::DREAMException(
                "Setting '%s': Expected value to be a scalar, but array was given.",
                name.c_str()
            );

        int dtype = PyArray_TYPE(ao);
        int_t v;
        if (dtype == NPY_INT) {
            v = reinterpret_cast<int*>(PyArray_DATA(ao))[0];
        } else if (dtype == NPY_LONG) {
            v = reinterpret_cast<long*>(PyArray_DATA(ao))[0];
        } else
            throw DREAM::DREAMException(
                "Setting '%s': Unrecognized data type of specified value: %d. Expected numpy integer array. (1)",
                name.c_str(), dtype
            );

        s->SetSetting(name, v);
    } else
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value. Expected integer.",
            name.c_str()
        );
}

/**
 * Load a setting as a real number from the
 * given Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_real(Settings *s, const string& name, PyObject *obj) {
    if (PyFloat_Check(obj)) {
        double d = PyFloat_AsDouble(obj);
        s->SetSetting(name, static_cast<real_t>(d));
    } else if (PyLong_Check(obj)) {
        long l = PyLong_AsLong(obj);
        s->SetSetting(name, static_cast<real_t>(l));
    } else if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        int ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        if (ndim != 1 || _dims[0] != 1)
            throw DREAM::DREAMException(
                "Setting '%s': Expected value to be a scalar, but array was given.",
                name.c_str()
            );

        int dtype = PyArray_TYPE(ao);
        real_t v;
        if (dtype == NPY_FLOAT) {
            v = reinterpret_cast<float*>(PyArray_DATA(ao))[0];
        } else if (dtype == NPY_DOUBLE) {
            v = reinterpret_cast<double*>(PyArray_DATA(ao))[0];
        } else if (dtype == NPY_INT) {
            v = reinterpret_cast<int*>(PyArray_DATA(ao))[0];
        } else if (dtype == NPY_LONG) {
            v = reinterpret_cast<long*>(PyArray_DATA(ao))[0];
        } else
            throw DREAM::DREAMException(
                "Setting '%s': Unrecognized data type of specified value: %d. Expected numpy integer array. (1)",
                name.c_str(), dtype
            );

        s->SetSetting(name, v);
    } else
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value. Expected real number.",
            name.c_str()
        );
}

/**
 * Load a setting as an integer array from
 * the given Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_int_array(Settings *s, const string& name, PyObject *obj) {
    if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        // Get dimensions of array
        int ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        if (ndim == 1 && _dims[0] == 1) {
        }

        // Convert to 'len_t' (needed for DREAM Settings API)
        len_t *dims = new len_t[ndim];
        for (int i = 0; i < ndim; i++)
            dims[i] = _dims[i];

        // Get array data type
        int dtype = PyArray_TYPE(ao);
        int_t *v;

        if (dtype == NPY_INT) {
            v = dreampy_convert<int_t,int>(
                reinterpret_cast<int*>(PyArray_DATA(ao)),
                ndim, _dims
            );
        } else if (dtype == NPY_LONG) {
            v = dreampy_convert<int_t,long>(
                reinterpret_cast<long*>(PyArray_DATA(ao)),
                ndim, _dims
            );
        } else
            throw DREAM::DREAMException(
                "Setting '%s': Unrecognized data type of specified value: %d. Expected numpy integer array. (1)",
                name.c_str(), dtype
            );

        // Set setting value
        s->SetSetting(name, ndim, dims, v);

        delete [] dims;
    } else if (PyList_Check(obj)) {
        Py_ssize_t n = PyList_Size(obj);

        int_t *v = new int_t[n];
        for (Py_ssize_t i = 0; i < n; i++) {
            PyObject *li = PyList_GetItem(obj, i);

            if (PyLong_Check(li))
                v[i] = PyLong_AsLong(li);
            else if (PyArray_IsAnyScalar(li)) {
                PyArray_Descr *type = PyArray_DescrFromType(NPY_INT64);
                PyArray_CastScalarToCtype(li, v+i, type);
                Py_DECREF(type);
            } else
                throw DREAM::DREAMException(
                    "Setting '%s': Unrecognized type of list element: %s",
                    name.c_str(), li->ob_type->tp_name
                );
        }

        len_t nel = n;
        s->SetSetting(name, 1, &nel, v);
    } else {
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value. Expected numpy integer array (2).",
            name.c_str()
        );
    }
}

/**
 * Load a setting as an real number array from
 * the given Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_real_array(Settings *s, const string& name, PyObject *obj) {
    if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        // Get dimensions of array
        int ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        // Convert to 'len_t' (needed for DREAM Settings API)
        len_t *dims = new len_t[ndim];
        for (int i = 0; i < ndim; i++)
            dims[i] = _dims[i];

        // Get array data type
        int dtype = PyArray_TYPE(ao);
        real_t *v;
        
        if (dtype == NPY_FLOAT) {
            v = dreampy_convert<real_t,float>(
                reinterpret_cast<float*>(PyArray_DATA(ao)),
                ndim, _dims
            );
        } else if (dtype == NPY_DOUBLE) {
            v = dreampy_convert<real_t,double>(
                reinterpret_cast<double*>(PyArray_DATA(ao)),
                ndim, _dims
            );
        } else if (dtype == NPY_LONG) {
            v = dreampy_convert<real_t,long>(
                reinterpret_cast<long*>(PyArray_DATA(ao)),
                ndim, _dims
            );
        } else
            throw DREAM::DREAMException(
                "Setting '%s': Unrecognized data type of specified value: %d. Expected numpy integer array. (1)",
                name.c_str(), dtype
            );

        // Set setting value
        s->SetSetting(name, ndim, dims, v);

        delete [] dims;
    } else if (PyList_Check(obj)) {
        Py_ssize_t n = PyList_Size(obj);

        real_t *v = new real_t[n];
        for (Py_ssize_t i = 0; i < n; i++) {
            PyObject *li = PyList_GetItem(obj, i);

            if (PyFloat_Check(li))
                v[i] = PyFloat_AsDouble(li);
            else
                throw DREAM::DREAMException(
                    "Setting '%s': Unrecognized type of list element: %s",
                    name.c_str(), li->ob_type->tp_name
                );
        }

        len_t nel = n;
        s->SetSetting(name, 1, &nel, v);
    } else {
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value. Expected numpy integer array (2).",
            name.c_str()
        );
    }
}

/**
 * Load a setting as a string from the given
 * Python object.
 *
 * s:    Settings object to assign value to.
 * name: Name of setting to assign.
 * obj:  Python object to load value from.
 */
void dreampy_load_string(Settings *s, const string& name, PyObject *obj) {
    if (PyUnicode_Check(obj)) {
        string str = PyUnicode_AsUTF8(obj);
        s->SetSetting(name, str);
    } else {
        throw DREAM::DREAMException(
            "Setting '%s': Unrecognized data type of specified value. Expected string. "
            "Data is of type '%s'.",
            name.c_str(), obj->ob_type->tp_name
        );
    }
}

