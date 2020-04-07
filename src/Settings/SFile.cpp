/**
 * LOAD SETTINGS VIA THE SOFTLIB SFILE INTERFACE
 * This module loads settings from a file using the softlib SFile interface,
 * i.e. from HDF5, MATLAB MAT or SDT files.
 */

#include <string>
#include <softlib/SFile.h>
#include "DREAM/IO.hpp"
#include "DREAM/Settings/SFile.hpp"
#include "DREAM/Settings/Settings.hpp"


using namespace DREAM;
using namespace std;


/**
 * Load settings from the given file.
 *
 * settings: Settings object to set.
 * filename: Name of file to load settings from.
 */
void DREAM::SettingsSFile::LoadSettings(Settings *settings, const string& filename) {
    SFile *sf = SFile::Create(filename, SFILE_MODE_READ);
    LoadSettings(settings, sf);
}
void DREAM::SettingsSFile::LoadSettings(Settings *settings, SFile *sf) {
    const map<string, Settings::setting_t*> allset = settings->GetSettings();
    vector<string> missing;

    // Load settings
    for (auto it = allset.begin(); it != allset.end(); it++) {
        if (sf->HasVariable(it->first)) {
            LoadSetting(it->first, it->second->type, it->second->ndims, sf, settings);
        } else if (it->second->mandatory) {
            missing.push_back(
                sf->filename + ": The mandatory setting '" + it->first + "' "
                "was not present in the file."
            );
        }
    }

    // Check if any mandatory settings are missing
    if (missing.size() > 0) {
        for (auto it = missing.begin(); it != missing.end(); it++)
            IO::PrintError(*it);

        throw SettingsException(
            "%s: Mandatory settings were not provided.",
            sf->filename.c_str()
        );
    }
}

/**
 * Load a single setting from the given SFile object
 * and store in the given Settings object.
 *
 * name: Name of setting to load.
 * sf:   SFile object to load setting with.
 * set:  Settings object to store setting in.
 */
void DREAM::SettingsSFile::LoadSetting(
    const string& name, enum Settings::setting_type type,
    const len_t ndims, SFile *sf, Settings *set
) {
    switch (type) {
        case Settings::SETTING_TYPE_BOOL: LoadBool(name, sf, set); break;
        case Settings::SETTING_TYPE_INT: LoadInteger(name, sf, set); break;
        case Settings::SETTING_TYPE_REAL: LoadReal(name, sf, set); break;
        case Settings::SETTING_TYPE_INT_ARRAY: LoadIntegerArray(name, ndims, sf, set); break;
        case Settings::SETTING_TYPE_REAL_ARRAY: LoadRealArray(name, ndims, sf, set); break;

        default:
            throw SettingsException(
                "%s: setting '%s': Unrecognized setting type: %d.",
                sf->filename.c_str(), name.c_str(), type
            );
    }
}

/**
 * Load a bool value from the given SFile into the
 * given Settings object.
 */
void DREAM::SettingsSFile::LoadBool(const string& name, SFile *sf, Settings *set) {
    int64_t v = sf->GetInt(name);
    set->SetSetting(name, (v != 0));
}

/**
 * Load integer value from the given SFile into the
 * given Settings object.
 */
void DREAM::SettingsSFile::LoadInteger(const string& name, SFile *sf, Settings *set) {
    int64_t v = sf->GetInt(name);
    set->SetSetting(name, (int_t)v);
}

/**
 * Load real value from the given SFile into the
 * given Settings object.
 */
void DREAM::SettingsSFile::LoadReal(const string& name, SFile *sf, Settings *set) {
    real_t v = (real_t)sf->GetScalar(name);
    set->SetSetting(name, v);
}

/**
 * Load array of integers from the given SFile into
 * the given Settings object.
 *
 * name:          Name of setting to load.
 * nExpectedDims: Expected number of dimensions in array.
 * sf:            SFile object to read from.
 * set:           Settings object to assign value to.
 */
void DREAM::SettingsSFile::LoadIntegerArray(
    const string& name, const len_t nExpectedDims,
    SFile *sf, Settings *set
) {
    int_t *v;
    len_t ndims;
    len_t *dims = new len_t[nExpectedDims];

    if (nExpectedDims != 1)
        throw SettingsException(
            "%s: setting '%s': Only 1D integer arrays are supported by the "
            "SFile interface at the moment.",
            sf->filename.c_str(), name.c_str()
        );

    // Load without having to convert data?
    if (typeid(int_t) == typeid(int64_t)) {
        sfilesize_t _ndims=nExpectedDims;
        sfilesize_t *_dims = new sfilesize_t[nExpectedDims];
        v = sf->GetIntList(name, _dims);

        ndims = _ndims;
        for (len_t i = 0; i < ndims; i++)
            dims[i] = (len_t)_dims[i];

        delete [] _dims;
    // int_t != int64_t  ==> convert data
    } else {
        sfilesize_t _ndims=nExpectedDims;
        sfilesize_t *_dims = new sfilesize_t[nExpectedDims];
        int64_t *d = sf->GetIntList(name, _dims);

        len_t ntot = 1;
        ndims = _ndims;
        for (len_t i = 0; i < ndims; i++) {
            dims[i] = (len_t)_dims[i];
            ntot *= dims[i];
        }

        v = new int_t[ntot];
        for (len_t i = 0; i < ntot; i++)
            v[i] = (real_t)d[i];

        delete [] _dims;
    }

    set->SetSetting(name, ndims, dims, v);

    delete [] dims;
}

/**
 * Load array of real values from the given SFile into
 * the given Settings object.
 *
 * name:          Name of setting to load.
 * nExpectedDims: Expected number of dimensions in array.
 * sf:            SFile object to read from.
 * set:           Settings object to assign value to.
 */
void DREAM::SettingsSFile::LoadRealArray(
    const string& name, const len_t nExpectedDims,
    SFile *sf, Settings *set
) {
    real_t *v;
    len_t ndims;
    len_t *dims = new len_t[nExpectedDims];

    // Load without having to convert data?
    if (typeid(real_t) == typeid(double)) {
        sfilesize_t _ndims;
        sfilesize_t *_dims = new sfilesize_t[nExpectedDims];
        v = sf->GetMultiArray_linear(name, nExpectedDims, _ndims, _dims);

        ndims = _ndims;
        for (len_t i = 0; i < ndims; i++)
            dims[i] = (len_t)_dims[i];

        delete [] _dims;
    } else {    // real_t != double  ==> convert data
        sfilesize_t _ndims;
        sfilesize_t *_dims = new sfilesize_t[nExpectedDims];
        double *d = sf->GetMultiArray_linear(name, nExpectedDims, _ndims, _dims);

        len_t ntot = 1;
        ndims = _ndims;
        for (len_t i = 0; i < ndims; i++) {
            dims[i] = (len_t)_dims[i];
            ntot *= dims[i];
        }

        v = new real_t[ntot];
        for (len_t i = 0; i < ntot; i++)
            v[i] = (real_t)d[i];

        delete [] _dims;
    }

    set->SetSetting(name, ndims, dims, v);

    delete [] dims;
}

