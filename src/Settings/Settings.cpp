/**
 * Implementation of the 'Settings' object.
 */

#include <algorithm>
#include <iostream>
#include <string>
#include "DREAM/Settings/Settings.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
Settings::Settings() {
}

/**
 * Destructor.
 */
Settings::~Settings() {
    for (auto it = settings.begin(); it != settings.end(); it++) {
        delete it->second;
    }
}

/******************************
 * INTERNAL METHODS           *
 ******************************/
/**
 * Define the named setting. This method first checks
 * that the setting has not been previously defined.
 *
 * name:         Name of setting to define.
 * desc:         Description of the setting.
 * defaultValue: Value assigned to setting by default.
 * type:         Data type of setting.
 * mandatory:    If 'true', the setting MUST be specified by
 *               the user.
 */
template<typename T>
void Settings::_DefineSetting(
    const string& name, const string& desc, const T& defaultValue,
    enum setting_type type, bool mandatory
) {
    if (settings.find(name) != settings.end())
        throw SettingsException("The setting '%s' has already been defined.", name.c_str());

    setting_t *s   = new setting_t;
    s->description = desc;
    s->type        = type;
    s->ndims       = 1;
    s->mandatory   = mandatory;

    T *a = new T;
    *a = defaultValue;
    s->value = a;

    settings[name] = s;
}

/**
 * Define the named setting. This method first checks
 * that the setting has not been previously defined.
 *
 * name:         Name of setting to define.
 * desc:         Description of the setting.
 * ndims:        Number of dimensions of array.
 * dims:         List of number of elements per dimension
 *               of array.
 * defaultValue: Value assigned to setting by default.
 * type:         Data type of setting.
 * mandatory:    If 'true', the setting MUST be specified by
 *               the user.
 */
template<typename T>
void Settings::_DefineSetting(
    const string& name, const string& desc,
    const len_t ndims, const len_t dims[], const T *defaultValue,
    enum setting_type type, bool mandatory
) {
    if (settings.find(name) != settings.end())
        throw SettingsException("The setting '%s' has already been defined.", name.c_str());

    setting_t *s   = new setting_t;
    s->description = desc;
    s->type        = type;
    s->ndims       = ndims;
    s->dims        = new len_t[ndims];
    s->mandatory   = mandatory;

    len_t ntot = 1;
    for (len_t i = 0; i < ndims; i++) {
        s->dims[i] = dims[i];
        ntot *= dims[i];
    }

    if (defaultValue != nullptr) {
        T *a = new T;
        for (len_t i = 0; i < ntot; i++)
            a[i] = defaultValue[i];
        
        s->value = a;
    } else s->value = nullptr;

    settings[name] = s;
}

/**
 * Remove the definition of the named setting from the settings object.
 * No check is made to verify that the setting has already been defined,
 * so if an attempt is made to remove a non-existant setting, an exception
 * will be thrown.
 *
 * name: Name of setting to remove.
 */
void Settings::UndefineSetting(const string& name) {
	settings.erase(name);
}

/**
 * Returns the setting with the given name, if defined.
 *
 * name:     Name of setting to retrieve.
 * type:     Setting data type.
 * markused: If 'true', marks the setting as "used".
 */
Settings::setting_t *Settings::_GetSetting(const string& name, enum setting_type type, bool markused) {
    auto it = settings.find(name);
    if (it == settings.end())
        throw SettingsException("The setting '%s' has not been defined.", name.c_str());

    setting_t *s = it->second;
    if (s->type != type)
        throw SettingsException(
            "The setting '%s' is not %s as expected. It is %s.",
            name.c_str(), GetTypeName(type), GetTypeName(s->type)
        );

    if (markused) s->used = true;

    return s;
}

/**
 * Returns the specified setting as an array of the
 * given type.
 *
 * name:          Name of setting load load.
 * nExpectedDims: Number of dimensions that we expect the
 *                array to have (i.e. 'ndims' contains
 *                this many elements)
 * ndims:         The number of elements in each array dimension.
 * markused:      If 'true', marks the setting as "used".
 */
template<typename T>
T *Settings::_GetArray(
    const string& name,
    const len_t nExpectedDims, len_t ndims[],
    enum setting_type type, bool markused
) {
    setting_t *s = _GetSetting(name, type, markused);

    if (nExpectedDims != s->ndims)
        throw SettingsException(
            "Setting '%s': Invalid number of dimensions of array. Expected "
            LEN_T_PRINTF_FMT " dimensions. Array has %d dimensions.",
            name.c_str(), nExpectedDims, s->ndims
        );

    for (len_t i = 0; i < nExpectedDims; i++)
        ndims[i] = s->dims[i];

    if (markused) s->used = true;

    return (T*)s->value;
}

/**
 * Set the value of a previously defined setting.
 *
 * name:  Name of setting to set.
 * value: Value to assign to setting.
 */
template<typename T>
void Settings::_SetSetting(
    const string& name, const T& value,
    enum setting_type type
) {
    auto it = settings.find(name);
    if (it == settings.end())
        throw SettingsException(
            "No setting named '%s' exists.", name.c_str()
        );

    if (it->second->type != type)
        throw SettingsException(
            "%s: The given value is %s, while %s was expected.",
            name.c_str(), GetTypeName(type), GetTypeName(it->second->type)
        );

    *((T*)(it->second->value)) = value;
}

/**
 * Set the value of a previously defined setting to
 * an array. NOTE: No new memory is allocated and
 * this Settings object assumes ownership of the array.
 * Hence, cleanup will be handled by this object.
 *
 * name:  Name of setting to set.
 * value: Value to assign to setting.
 */
template<typename T>
void Settings::_SetSetting(
    const string& name, len_t ndims,
    const len_t dims[], T *value,
    enum setting_type type
) {
    auto it = settings.find(name);
    if (it == settings.end())
        throw SettingsException(
            "No setting named '%s' exists.", name.c_str()
        );

    if (it->second->type != type)
        throw SettingsException(
            "Setting '%s': The given value is %s, while %s was expected.",
            name.c_str(), GetTypeName(type), GetTypeName(it->second->type)
        );
    
    if (it->second->ndims != ndims)
        throw SettingsException(
            "Setting: '%s': The given value has an invalid dimensionality: "
            LEN_T_PRINTF_FMT ". Expected dimensionality: " LEN_T_PRINTF_FMT,
            name.c_str(), ndims, it->second->ndims
        );

    for (len_t i = 0; i < ndims; i++)
        it->second->dims[i] = dims[i];

    it->second->value = value;
}

/***************************
 * PUBLIC METHODS          *
 ***************************/
void Settings::DefineSetting(const string& name, const string& desc, bool defaultValue, bool mandatory)
{ this->_DefineSetting<bool>(name, desc, defaultValue, SETTING_TYPE_BOOL, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, int_t defaultValue, bool mandatory)
{ this->_DefineSetting<int_t>(name, desc, defaultValue, SETTING_TYPE_INT, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, real_t defaultValue, bool mandatory)
{ this->_DefineSetting<real_t>(name, desc, defaultValue, SETTING_TYPE_REAL, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, const string& defaultValue, bool mandatory)
{ this->_DefineSetting<string>(name, desc, defaultValue, SETTING_TYPE_STRING, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, len_t n, const int_t *defaultValue, bool mandatory)
{ this->_DefineSetting<int_t>(name, desc, 1, &n, defaultValue, SETTING_TYPE_INT_ARRAY, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, len_t n, const len_t dims[], const int_t *defaultValue, bool mandatory)
{ this->_DefineSetting<int_t>(name, desc, n, dims, defaultValue, SETTING_TYPE_INT_ARRAY, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, len_t n, const real_t *defaultValue, bool mandatory)
{ this->_DefineSetting<real_t>(name, desc, 1, &n, defaultValue, SETTING_TYPE_REAL_ARRAY, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, len_t n, const len_t dims[], const real_t *defaultValue, bool mandatory)
{ this->_DefineSetting<real_t>(name, desc, n, dims, defaultValue, SETTING_TYPE_REAL_ARRAY, mandatory); }

void Settings::DefineSetting(const string& name, const string& desc, void *defaultValue, bool mandatory)
{ this->_DefineSetting<void*>(name, desc, defaultValue, SETTING_TYPE_ADDRESS, mandatory); }

/**
 * Returns 'true' if the named setting has been defined.
 */
bool Settings::HasSetting(const string& name) {
    auto it = settings.find(name);
    return (it != settings.end());
}

/**
 * Returns the data type of the specified setting.
 *
 * name: Name of setting to return type of.
 */
enum Settings::setting_type Settings::GetType(const string& name) {
    auto it = settings.find(name);
    if (it == settings.end())
        throw SettingsException("The setting '%s' has not been defined.", name.c_str());

    setting_t *s = it->second;
    return s->type;
}

/**
 * Returns the specified setting as an address.
 */
void *Settings::GetAddress(const string& name, bool markused) {
    return *((void**)(_GetSetting(name, SETTING_TYPE_ADDRESS, markused)->value));
}

/**
 * Returns the specified setting as a bool.
 */
bool Settings::GetBool(const string& name, bool markused) {
    return *((bool*)(_GetSetting(name, SETTING_TYPE_BOOL, markused)->value));
}

/**
 * Returns the specified setting as an integer.
 */
int_t Settings::GetInteger(const string& name, bool markused) {
    return *((int_t*)(_GetSetting(name, SETTING_TYPE_INT, markused)->value));
}

/**
 * Returns the specified setting as a real number.
 */
real_t Settings::GetReal(const string& name, bool markused) {
    return *((real_t*)(_GetSetting(name, SETTING_TYPE_REAL, markused)->value));
}

/**
 * Returns the specified setting as a string.
 */
const string Settings::GetString(const string& name, bool markused) {
    return *((string*)(_GetSetting(name, SETTING_TYPE_STRING, markused)->value));
}

/**
 * Returns the specified setting as a list of strings. The list is
 * stored as a single string in this Settings object, and so this
 * routine retrieves the Setting as a regular string first and the
 * splits it based on the specified delimiter.
 *
 * name:     Name of setting to retrieve.
 * delim:    Delimiter character used to separate list elements.
 * markused: If 'true', marks the setting as used.
 */
vector<string> Settings::GetStringList(
    const string& name, const char delim, bool markused
) {
    const string s = GetString(name, markused);
    vector<string> list;
    len_t si = 0;
    const len_t sl = s.size();

    if (sl > 0)
        list.push_back("");

    for (len_t i = 0; i < sl; i++) {
        if (s[i] == delim) {
            si++;
            if (i+1 < sl)
                list.push_back("");
        } else
            list[si] += s[i];
    }

    return list;
}

/**
 * Returns the specified setting as an integer array.
 *
 * name:          Name of setting load load.
 * nExpectedDims: Number of dimensions that we expect the
 *                array to have (i.e. 'ndims' contains
 *                this many elements)
 * ndims:         The number of elements in each array dimension.
 * markused:      If 'true', marks the settings as "used".
 */
const int_t *Settings::GetIntegerArray(
    const string& name, const len_t nExpectedDims, len_t ndims[], bool markused
) {
    return _GetArray<int_t>(name, nExpectedDims, ndims, SETTING_TYPE_INT_ARRAY, markused);
}

/**
 * Returns the specified setting as a real array.
 *
 * name:          Name of setting load load.
 * nExpectedDims: Number of dimensions that we expect the
 *                array to have (i.e. 'ndims' contains
 *                this many elements)
 * ndims:         The number of elements in each array dimension.
 */
const real_t *Settings::GetRealArray(
    const string& name, const len_t nExpectedDims, len_t ndims[], bool markused
) {
    return _GetArray<real_t>(name, nExpectedDims, ndims, SETTING_TYPE_REAL_ARRAY, markused);
}

/**
 * Mark the specified setting as "used".
 *
 * name: Name of setting to mark as used.
 */
void Settings::MarkUsed(const std::string& name) {
    auto it = settings.find(name);
    if (it == settings.end())
        throw SettingsException("The setting '%s' has not been defined.", name.c_str());
    else
        it->second->used = true;
}


void Settings::SetSetting(const string& name, bool value)
{ this->_SetSetting(name, value, SETTING_TYPE_BOOL); }

void Settings::SetSetting(const string& name, int_t value)
{ this->_SetSetting(name, value, SETTING_TYPE_INT); }

void Settings::SetSetting(const string& name, real_t value)
{ this->_SetSetting(name, value, SETTING_TYPE_REAL); }

void Settings::SetSetting(const string& name, const string& value)
{ this->_SetSetting(name, value, SETTING_TYPE_STRING); }

void Settings::SetSetting(const string& name, const len_t n, int_t *value)
{ this->_SetSetting(name, 1, &n, value, SETTING_TYPE_INT_ARRAY); }

void Settings::SetSetting(const string& name, const len_t ndims, const len_t dims[], int_t *value)
{ this->_SetSetting(name, ndims, dims, value, SETTING_TYPE_INT_ARRAY); }

void Settings::SetSetting(const string& name, const len_t n, real_t *value)
{ this->_SetSetting(name, 1, &n, value, SETTING_TYPE_REAL_ARRAY); }

void Settings::SetSetting(const string& name, const len_t ndims, const len_t dims[], real_t *value)
{ this->_SetSetting(name, ndims, dims, value, SETTING_TYPE_REAL_ARRAY); }

void Settings::SetSetting(const string& name, void *value)
{ this->_SetSetting(name, value, SETTING_TYPE_ADDRESS); }

/**
 * Print a list of all available settings.
 */
void Settings::DisplaySettings() {
    for (auto it = settings.begin(); it != settings.end(); it++) {
        setting_t *s = it->second;
        printf("%-40s -- %s\n", it->first.c_str(), s->description.c_str());
    }
}

