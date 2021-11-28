/**
 * LOAD/SAVE SETTINGS VIA THE SOFTLIB SFILE INTERFACE
 * This module loads settings from a file using the softlib SFile interface,
 * i.e. from HDF5, MATLAB MAT or SDT files.
 */

#include <map>
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
            IO::PrintError("%s", it->c_str());

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
        case Settings::SETTING_TYPE_STRING: LoadString(name, sf, set); break;
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
 * Load string from the given SFile into the given
 * Settings object.
 */
void DREAM::SettingsSFile::LoadString(const string& name, SFile *sf, Settings *set) {
    string v = sf->GetString(name);
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

/**
 * Creates any groups in the 'name' which do not already exist.
 *
 * name:   Name of setting to be saved (including full path relative to 'path').
 * path:   Path relative to root of file to save setting in (assumed to exist).
 * groups: List of groups which have already been created.
 */
void DREAM::SettingsSFile::CreateGroup(
    const string& name, const string& path, vector<string> &groups,
    SFile *sf
) {
    size_t slash = name.find('/'), start = 0;
    string group = "";
    while (slash != string::npos) {
        group += name.substr(start, slash-start);

        // If the group does not exist in 'groups'...
        if (std::find(groups.begin(), groups.end(), group) == std::end(groups)) {
            sf->CreateStruct(path+group);
            groups.push_back(group);
        }

        group += '/';
        start = slash+1;
        slash = name.find('/', start);
    }
}

/**
 * Save the given Settings object using the SFile interface.
 *
 * settings: Settings object to write.
 * sf:       SFile object to write to.
 * path:     Path in output file to write settings to.
 */
void DREAM::SettingsSFile::SaveSettings(
    Settings *s, SFile *sf, const string& path
) {
    string group = path;
    if (path.back() != '/')
        group += '/';

    const map<string, Settings::setting_t*> smap = s->GetSettings();
    vector<string> groups;

    for (auto const& [name, set] : smap) {
        // Ignore unused settings
        if (not set->used)
            continue;

        string sname;
        if (name[0] == '/')
            sname = name.substr(1);
        else
            sname = name;

        // Construct full name of setting
        string fullname = group + sname;

        CreateGroup(sname, group, groups, sf);

        // Should we create a new group?
        /*auto slash = name.find('/');
        if (slash != string::npos) {
            string groupname = name.substr(0, slash);

            if (std::find(groups.begin(), groups.end(), groupname) == 
        }*/

        switch (set->type) {
            case Settings::SETTING_TYPE_BOOL:
                SettingsSFile::SaveBool(fullname, s->GetBool(name), sf);
                break;
            case Settings::SETTING_TYPE_INT:
                SettingsSFile::SaveInteger(fullname, s->GetInteger(name), sf);
                break;
            case Settings::SETTING_TYPE_INT_ARRAY: {
                len_t *dims = new len_t[set->ndims];
                const int_t *arr = s->GetIntegerArray(name, set->ndims, dims);

                SettingsSFile::SaveIntegerArray(fullname, arr, set->ndims, dims, sf);

                delete [] dims;
            } break;
            case Settings::SETTING_TYPE_REAL:
                SettingsSFile::SaveReal(fullname, s->GetReal(name), sf);
                break;
            case Settings::SETTING_TYPE_REAL_ARRAY: {
                len_t *dims = new len_t[set->ndims];
                const real_t *arr = s->GetRealArray(name, set->ndims, dims);

                SettingsSFile::SaveRealArray(fullname, arr, set->ndims, dims, sf);

                delete [] dims;
            } break;
            case Settings::SETTING_TYPE_STRING:
                SettingsSFile::SaveString(fullname, s->GetString(name), sf);
                break;
            
            default:
                throw SettingsException(
                    "SettingsSFile: Unrecognized setting type: " LEN_T_PRINTF_FMT,
                    set->type
                );
        }
    }
}

/**
 * Save a value of type 'bool' to the given SFile object.
 *
 * name: Full path to setting in SFile object.
 * b:    Value to write.
 * sf:   SFile object to write to.
 */
void DREAM::SettingsSFile::SaveBool(
    const string& name, bool b, SFile *sf
) {
    int64_t v = b?1:0;
    sf->WriteInt64List(name, &v, 1);
}

/**
 * Save a value of type 'real' to the given SFile object.
 *
 * name: Full path to setting in SFile object.
 * r:    Value to write.
 * sf:   SFile object to write to.
 */
void DREAM::SettingsSFile::SaveReal(
    const string& name, real_t r, SFile *sf
) {
    double v = r;
    sf->WriteScalar(name, v);
}

/**
 * Save a value of type 'int_t' to the given SFile object.
 *
 * name: Full path to setting in SFile object.
 * i:    Value to write.
 * sf:   SFile object to write to.
 */
void DREAM::SettingsSFile::SaveInteger(
    const string& name, int_t i, SFile *sf
) {
    int64_t v = i;
    sf->WriteInt64List(name, &v, 1);
}

/**
 * Save a value of type 'real_t*' to the given SFile object.
 *
 * name:  Full path to setting in SFile object.
 * arr:   Array to write to SFile.
 * ndims: Number of dimensions of array.
 * dims:  Size of each dimension of array.
 * sf:    SFile object to write to.
 */
void DREAM::SettingsSFile::SaveRealArray(
    const string& name, const real_t *arr,
    const len_t ndims, const len_t *dims,
    SFile *sf
) {
    sfilesize_t *_dims = new sfilesize_t[ndims];
    for (len_t i = 0; i < ndims; i++)
        _dims[i] = dims[i];

    sf->WriteMultiArray(name, arr, ndims, _dims);

    delete [] _dims;
}

/**
 * Save a value of type 'int_t*' to the given SFile object.
 *
 * name:  Full path to setting in SFile object.
 * arr:   Array to write to SFile.
 * ndims: Number of dimensions of array.
 * dims:  Size of each dimension of array.
 * sf:    SFile object to write to.
 */
void DREAM::SettingsSFile::SaveIntegerArray(
    const string& name, const int_t *arr,
    const len_t ndims, const len_t *dims,
    SFile *sf
) {
    sfilesize_t *_dims = new sfilesize_t[ndims];
    for (len_t i = 0; i < ndims; i++)
        _dims[i] = dims[i];

    //sf->WriteMultiInt64Array(name, arr, ndims, _dims);

    delete [] _dims;
}

/**
 * Save a value of type 'string' to the given SFile object.
 *
 * name: Full path to the setting in the SFile object.
 * str:  String to save.
 * sf:   SFile object to save value to.
 */
void DREAM::SettingsSFile::SaveString(
    const string& name, const string& str, SFile *sf
) {
    sf->WriteString(name, str);
}

