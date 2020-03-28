#ifndef _DREAM_SETTINGS_HPP
#define _DREAM_SETTINGS_HPP

#include <map>
#include <string>
#include "FVM/config.h"
#include "FVM/FVMException.hpp"

namespace DREAM {
    class SettingsException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        SettingsException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Settings");
        }
    };

    class Settings {
    public:
        enum setting_type {
            SETTING_TYPE_BOOL,
            SETTING_TYPE_INT,
            SETTING_TYPE_INT_ARRAY,
            SETTING_TYPE_REAL,
            SETTING_TYPE_REAL_ARRAY
        };
        typedef struct _setting {
            std::string description;
            enum setting_type type;
            len_t ndims, *dims=nullptr;
            void *value;

            ~_setting() {
                switch (this->type) {
                    case SETTING_TYPE_BOOL: delete (bool*)value; break;
                    case SETTING_TYPE_INT:  delete (int*)value; break;
                    case SETTING_TYPE_REAL: delete (real_t*)value; break;
                    case SETTING_TYPE_INT_ARRAY: delete [] (int*)value; break;
                    case SETTING_TYPE_REAL_ARRAY: delete [] (real_t*)value; break;

                    default: break;
                }

                if (dims != nullptr)
                    delete [] dims;
            }
        } setting_t;

    private:
        std::map<std::string, setting_t*> settings;

        const char *GetTypeName(enum setting_type s) {
            switch (s) {
                case SETTING_TYPE_BOOL: return "a boolean";
                case SETTING_TYPE_INT: return "an integer";
                case SETTING_TYPE_INT_ARRAY: return "an array of integers";
                case SETTING_TYPE_REAL: return "a real number";
                case SETTING_TYPE_REAL_ARRAY: return "an array of real numbers";

                default: throw SettingsException("Unrecognized setting type: %d.", s);
            }
        }

        template<typename T>
        void _DefineSetting(const std::string& name, const std::string& desc, T& defaultValue, enum setting_type type);

        template<typename T>
        void _DefineSetting(const std::string& name, const std::string& desc, const len_t ndims, const len_t dims[], const T *defaultValue, enum setting_type type);

        template<typename T>
        T *_GetArray(const std::string&, const len_t, const len_t[], enum setting_type);
        setting_t *_GetSetting(const std::string&, enum setting_type);

        template<typename T>
        void _SetSetting(const std::string& name, const T& value, enum setting_type type);

        template<typename T>
        void _SetSetting(const std::string& name, const len_t ndims, const len_t dims[], T *value, enum setting_type type);

    public:
        Settings();
        ~Settings();

        void DefineSetting(const std::string& name, const std::string& desc, bool defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, int defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, len_t n, const int *defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, len_t ndims, const len_t dims[], const int *defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, real_t defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, len_t n, const real_t *defaultValue);
        void DefineSetting(const std::string& name, const std::string& desc, len_t ndims, const len_t dims[], const real_t *defaultValue);

        // GETTERS
        bool GetBool(const std::string&);
        int GetInteger(const std::string&);
        int *GetIntegerArray(const std::string& name, const len_t nExpectedDims, const len_t ndims[]);
        real_t GetReal(const std::string&);
        real_t *GetRealArray(const std::string& name, const len_t nExpectedDims, const len_t ndims[]);

        // SETTERS
        void SetSetting(const std::string& name, bool value);
        void SetSetting(const std::string& name, int value);
        void SetSetting(const std::string& name, len_t n, int *value);
        void SetSetting(const std::string& name, len_t ndims, const len_t dims[], int *value);
        void SetSetting(const std::string& name, real_t value);
        void SetSetting(const std::string& name, len_t n, real_t *value);
        void SetSetting(const std::string& name, len_t ndims, const len_t dims[], real_t *value);
    };
}

#endif/*_DREAM_SETTINGS_HPP*/
