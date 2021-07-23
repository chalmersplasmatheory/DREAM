/**
 * Implementation of a softlib SFile interface to Python dictionaries.
 */

#include <string>
#include "pyface/numpy.h"
#include "pyface/SFile_Python.hpp"
#include "DREAM/DREAMException.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "pyface/dreamtypes.hpp"
#include "pyface/settings.hpp"


using namespace std;


SFile_Python::SFile_Python() {
    this->dict = PyDict_New();
}
SFile_Python::~SFile_Python() {}

void SFile_Python::Close() { }


/**
 * Returns the dict with the given name from the given
 * Python dict.
 *
 * obj:  Container dictionary to take the dict from.
 * name: Name of dict to fetch.
 */
PyObject *SFile_Python::getDict(PyObject *obj, const string& name) {
    PyObject *o = PyDict_GetItemString(obj, name.c_str());

    if (!PyDict_Check(o))
        throw DREAM::DREAMException(
            "Key '%s': Object is not a dict.",
            name.c_str()
        );

    return o;
}

/**
 * Returns the closest parent 'dict' Python object of the
 * named key (which includes a full path from the root element).
 */
PyObject *SFile_Python::getParentStruct(const string& name, string& dsetname) {
    string s(name);
    size_t i = s.find('/');
    PyObject *obj = this->dict;

    while (i != string::npos) {
        // Ignore if '/' is the leading character...
        if (i > 0) {
            string d = s.substr(0, i);
            obj = getDict(obj, d);
        }

        s = s.substr(i+1);
        i = s.find('/');
    }

    dsetname = s;

    return obj;
}

/**
 * Returns the object corresponding to that named by the full
 * object path 'name'.
 */
PyObject *SFile_Python::getObjectByName(const std::string& name) {
    string dsetname;
    PyObject *dict = getParentStruct(name, dsetname);
    return PyDict_GetItemString(dict, dsetname.c_str());
}

/**
 * Create a new empty dictionary and insert in the specified
 * location of the output dict.
 *
 * name: Full path to the new dictionary.
 */
void SFile_Python::CreateStruct(const string& name) {
    string dname;
    PyObject *par = getParentStruct(name, dname);

    PyObject *newdict = PyDict_New();
    PyDict_SetItemString(par, dname.c_str(), newdict);
}

/**
 * Load attribute value from the named dataset (as a floating-point value).
 *
 * datasetname: Name of dataset to which the attribute applies.
 * name:        Name of attribute to load.
 */
double SFile_Python::GetAttributeScalar(const string&, const string&) {
    // TODO
    throw DREAM::NotImplementedException(
        "SFile_Python::GetAttributeScalar() has not been implemented yet."
    );
}

/**
 * Load attribute value from the named dataset (as a string).
 *
 * datasetname: Name of dataset to which the attribute applies.
 * name:        Name of attribute to load.
 */
string SFile_Python::GetAttributeString(const string&, const string&) {
    // TODO
    throw DREAM::NotImplementedException(
        "SFile_Python::GetAttributeString() has not been implemented yet."
    );
}

/**
 * Returns the data type of the named dataset.
 *
 * name: Name of dataset to check the type for.
 */
enum SFile::sfile_data_type SFile_Python::GetDataType(const string& name, enum SFile::sfile_data_type) {
    PyObject *obj = getObjectByName(name);

    if (PyFloat_Check(obj))
        return SFile::SFILE_DATA_DOUBLE;
    else if (PyLong_Check(obj))
        return SFile::SFILE_DATA_INT64;
    else if (PyUnicode_Check(obj))
        return SFile::SFILE_DATA_STRING;
    else
        return SFile::SFILE_DATA_UNDEFINED;
}

/**
 * Returns true if the named dataset/dict exists in the current
 * Python dictionary.
 */
bool SFile_Python::HasVariable(const string& name) {
    try {
        getObjectByName(name);
        return true;
    } catch (DREAM::DREAMException& ex) {
        return false;
    }
}

/**
 * Returns a multidimensional array as a linear vector.
 *
 * name:   Name of dataset to load.
 * nndims: Maximum allowed number of dimensions.
 * ndims:  Number of actual dimensions of array (return).
 * dims:   Size of each dimension on the array (return).
 */
template<typename T>
T *SFile_Python::getMultiArrayLinear(
    const string& name, const sfilesize_t nndims,
    sfilesize_t& ndims, sfilesize_t *dims
) {
    PyObject *obj = getObjectByName(name);

    if (PyArray_Check(obj)) {
        PyArrayObject *ao = reinterpret_cast<PyArrayObject*>(obj);

        sfilesize_t _ndim = PyArray_NDIM(ao);
        npy_intp *_dims = PyArray_DIMS(ao);

        if (_ndim > nndims)
            return nullptr;

        // Store numpy '_dims' array in returned 'dims'
        ndims = _ndim;
        for (sfilesize_t i = 0; i < ndims; i++)
            dims[i] = _dims[i];

        int dtype = PyArray_TYPE(ao);
        T *v;

        if (dtype == NPY_FLOAT) {
            v = dreampy_convert<T, float>(
                reinterpret_cast<float*>(PyArray_DATA(ao)),
                ndims, _dims
            );
        } else if (dtype == NPY_DOUBLE) {
            v = dreampy_convert<T, double>(
                reinterpret_cast<double*>(PyArray_DATA(ao)),
                ndims, _dims
            );
        } else if (dtype == NPY_LONG) {
            v = dreampy_convert<T, long>(
                reinterpret_cast<long*>(PyArray_DATA(ao)),
                ndims, _dims
            );
        } else
            throw DREAM::DREAMException(
                "Key '%s': Unrecognized data type of specified value: %d. Expected numpy integer array (1).",
                name.c_str(), dtype
            );

        return v;
    } else if (PyList_Check(obj)) {
        Py_ssize_t n = PyList_Size(obj);

        ndims = 1;
        dims[0] = n;

        T *v = new T[n];
        for (Py_ssize_t i = 0; i < n; i++) {
            PyObject *li = PyList_GetItem(obj, i);

            if (PyFloat_Check(li))
                v[i] = PyFloat_AsDouble(li);
            else
                throw DREAM::DREAMException(
                    "Key '%s': Unrecognized data type of specified value. Expected numpy integer array (2).",
                    name.c_str()
                );
        }

        return v;
    } else
        throw DREAM::DREAMException(
            "Key '%s': Unrecognized data type of specified value. Expected numpy integer array (3).",
            name.c_str()
        );
}

double *SFile_Python::GetMultiArray_linear(
    const string& name, const sfilesize_t nndims,
    sfilesize_t& ndims, sfilesize_t *dims
) {
    return getMultiArrayLinear<double>(name, nndims, ndims, dims);
}

/**
 * Read the named dataset as a string.
 */
string SFile_Python::GetString(const string& name) {
    PyObject *obj = getObjectByName(name);

    if (!PyUnicode_Check(obj))
        throw DREAM::DREAMException(
            "Key '%s': The data type is not 'str' as expected.",
            name.c_str()
        );

    string s = PyUnicode_AsUTF8(obj);

    return s;
}

/**
 * Dummy implementation for 'Open()'. Throw an exception since
 * this method should not be used for this SFile type.
 */
void SFile_Python::Open(const string&, enum sfile_mode) {
    throw DREAM::DREAMException(
        "SFile_Python::Open() cannot be used."
    );
}

/**
 * Write a dataset attribute consisting of a scalar floating-point
 * value to the named dataset.
 *
 * dataset: Name of dataset to set attribute for.
 * name:    Name of attribute.
 * v:       Scalar value to assign.
 */
void SFile_Python::WriteAttribute_scalar(
    const string& dataset, const string& name, const double v
) {
    string att = dataset + "@@";

    // Attributes are stored in dict's named the same as the
    // dataset, but with '@@' appended to it.
    PyObject *par=nullptr;
    try {
        // Try to open if attribute dict already exists...
        //par = getParentStruct(att+"/"+name, ts);
        par = getObjectByName(att);
    } catch (DREAM::DREAMException &ex) {}

    if (par == nullptr) {
        // Otherwise, first create it...
        CreateStruct(att);
        //par = getParentStruct(att+"/"+name, ts);
        par = getObjectByName(att);
    }

    PyDict_SetItemString(par, name.c_str(), PyFloat_FromDouble(v));
}
void SFile_Python::WriteAttribute_string(
    const string& dataset, const string& name, const string& v
) {
    string att = dataset + "@@";

    // Attributes are stored in dict's named the same as the
    // dataset, but with '@@' appended to it.
    PyObject *par = nullptr;
    try {
        // Try to open if attribute dict already exists...
        par = getObjectByName(att);
    } catch (DREAM::DREAMException &ex) {}

    if (par == nullptr) {
        // Otherwise, first create it...
        CreateStruct(att);
        //par = getParentStruct(att+"/"+name, ts);
        par = getObjectByName(att);
    }

    PyObject *s = PyUnicode_FromString(v.c_str());
    PyDict_SetItemString(par, name.c_str(), s);
}

/**
 * Write an array representing an image to the output.
 *
 * name: Name of dataset to store image in.
 * img:  Image data.
 * n:    Number of pixels in each dimension (same in both image dimensions).
 */
void SFile_Python::WriteImage(const string& name, const double *const* img, sfilesize_t n) {
    WriteArray(name, img, n, n);
}

/**
 * Write a multi-dimensional array to the output.
 *
 * name:  Name of dataset to store data in.
 * data:  Array data.
 * ndims: Number of dimensions of array.
 * dims:  List of size of each dimension of the array.
 * dtype: NumPy type of data.
 */
template<typename T>
void SFile_Python::writeArray(
    const string& name, const T* data, const sfilesize_t ndims,
    const sfilesize_t *dims, const int dtype
) {
    npy_intp *_dims = new npy_intp[ndims];
    sfilesize_t nel = 1;
    for (sfilesize_t i = 0; i < ndims; i++) {
        _dims[i] = dims[i];
        nel *= dims[i];
    }

    PyObject *arr = PyArray_SimpleNew(ndims, _dims, dtype);
    // Insert data...
    T *p = reinterpret_cast<T*>(
        PyArray_DATA(reinterpret_cast<PyArrayObject*>(arr))
    );
    for (sfilesize_t i = 0; i < nel; i++)
        p[i] = data[i];

    // Insert into dictionary...
    string dname;
    PyObject *par = getParentStruct(name, dname);
    PyDict_SetItemString(par, dname.c_str(), arr);
}

void SFile_Python::WriteMultiArray(
    const string& name, const double* data, const sfilesize_t ndims,
    const sfilesize_t *dims
) {
    writeArray<double>(name, data, ndims, dims, NPY_DOUBLE);
}

/**
 * Write a string to the output.
 *
 * name: Name of dataset to store string in.
 * str:  String to insert.
 */
void SFile_Python::WriteString(
    const string& name, const string& str
) {
    string dname;
    PyObject *par = getParentStruct(name, dname);
    // Construct Python string object...
    PyObject *s = PyUnicode_FromString(str.c_str());

    PyDict_SetItemString(par, dname.c_str(), s);
}

/**
 * Get a 1D array of the specified type from the dict.
 */
template<typename T>
T *SFile_Python::get1DArray(const string& name, sfilesize_t *size) {
    sfilesize_t adims;
    return getMultiArrayLinear<T>(name, 1, adims, size);
}

double *SFile_Python::GetDoubles1D(const string& name, sfilesize_t *size) {
    return get1DArray<double>(name, size);
}

int32_t *SFile_Python::GetInt32_1D(const string& name, sfilesize_t *size) {
    return get1DArray<int32_t>(name, size);
}

int64_t *SFile_Python::GetInt64_1D(const string& name, sfilesize_t *size) {
    return get1DArray<int64_t>(name, size);
}

uint32_t *SFile_Python::GetUInt32_1D(const string& name, sfilesize_t *size) {
    return get1DArray<uint32_t>(name, size);
}

uint64_t *SFile_Python::GetUInt64_1D(const string& name, sfilesize_t *size) {
    return get1DArray<uint64_t>(name, size);
}


/**
 * Get a 2D array of the specified type from the dict.
 */
template<typename T>
T **SFile_Python::get2DArray(const string& name, sfilesize_t *dims) {
    sfilesize_t adims;
    T *arr = getMultiArrayLinear<T>(name, 2, adims, dims);

    if (adims != 2)
        throw DREAM::DREAMException(
            "'%s': Invalid dimensions of array: %llu. Expected 2.",
            name.c_str(), adims
        );

    T **oarr = new T*[dims[0]];
    oarr[0] = arr;
    for (sfilesize_t i = 1; i < dims[0]; i++)
        oarr[i] = oarr[i-1] + dims[1];

    return oarr;
}

/**
 * Get a 2D array of doubles from the dict.
 */
double **SFile_Python::GetDoubles(const string& name, sfilesize_t *dims) {
    return get2DArray<double>(name, dims);
}

int32_t **SFile_Python::GetInt32_2D(const string& name, sfilesize_t *dims) {
    return get2DArray<int32_t>(name, dims);
}

int64_t **SFile_Python::GetInt64_2D(const string& name, sfilesize_t *dims) {
    return get2DArray<int64_t>(name, dims);
}

uint32_t **SFile_Python::GetUInt32_2D(const string& name, sfilesize_t *dims) {
    return get2DArray<uint32_t>(name, dims);
}

uint64_t **SFile_Python::GetUInt64_2D(const string& name, sfilesize_t *dims) {
    return get2DArray<uint64_t>(name, dims);
}

/**
 * Write a 2D array to the dictionary.
 */
void SFile_Python::WriteArray(
    const string& name, const double *const* data,
    sfilesize_t rows, sfilesize_t cols
) {
    sfilesize_t dims[2] = {rows, cols};
    return writeArray<double>(name, data[0], 2, dims, NPY_DOUBLE);
}

void SFile_Python::WriteInt32Array(
    const string& name, const int32_t *const* data,
    sfilesize_t cols, sfilesize_t rows
) {
    sfilesize_t dims[2] = {rows, cols};
    return writeArray<int32_t>(name, data[0], 2, dims, NPY_INT32);
}

void SFile_Python::WriteInt64Array(
    const string& name, const int64_t *const* data,
    sfilesize_t cols, sfilesize_t rows
) {
    sfilesize_t dims[2] = {rows, cols};
    return writeArray<int64_t>(name, data[0], 2, dims, NPY_INT64);
}

void SFile_Python::WriteUInt32Array(
    const string& name, const uint32_t *const* data,
    sfilesize_t cols, sfilesize_t rows
) {
    sfilesize_t dims[2] = {rows, cols};
    return writeArray<uint32_t>(name, data[0], 2, dims, NPY_UINT32);
}

void SFile_Python::WriteUInt64Array(
    const string& name, const uint64_t *const* data,
    sfilesize_t cols, sfilesize_t rows
) {
    sfilesize_t dims[2] = {rows, cols};
    return writeArray<uint64_t>(name, data[0], 2, dims, NPY_UINT64);
}

void SFile_Python::WriteList(
    const string& name, const double *data, sfilesize_t size
) {
    writeArray<double>(name, data, 1, &size, NPY_DOUBLE);
}

void SFile_Python::WriteInt32List(
    const string& name, const int32_t *data, sfilesize_t size
) {
    writeArray<int32_t>(name, data, 1, &size, NPY_INT32);
}

void SFile_Python::WriteInt64List(
    const string& name, const int64_t *data, sfilesize_t size
) {
    writeArray<int64_t>(name, data, 1, &size, NPY_INT64);
}

void SFile_Python::WriteUInt32List(
    const string& name, const uint32_t *data, sfilesize_t size
) {
    writeArray<uint32_t>(name, data, 1, &size, NPY_UINT32);
}

void SFile_Python::WriteUInt64List(
    const string& name, const uint64_t *data, sfilesize_t size
) {
    writeArray<uint64_t>(name, data, 1, &size, NPY_UINT64);
}

