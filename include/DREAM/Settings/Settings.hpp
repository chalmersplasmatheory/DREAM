#ifndef _DREAM_SETTINGS_HPP
#define _DREAM_SETTINGS_HPP

#include <map>
#include <string>
#include <vector>
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
            SETTING_TYPE_REAL_ARRAY,
            SETTING_TYPE_STRING
        };
        typedef struct _setting {
            std::string description;
            enum setting_type type;
            len_t ndims, *dims=nullptr;
            void *value;
            bool mandatory=false;
            bool used=false;

            ~_setting() {
                switch (this->type) {
                    case SETTING_TYPE_BOOL: delete (bool*)value; break;
                    case SETTING_TYPE_INT:  delete (int_t*)value; break;
                    case SETTING_TYPE_REAL: delete (real_t*)value; break;
                    case SETTING_TYPE_INT_ARRAY: delete [] (int_t*)value; break;
                    case SETTING_TYPE_REAL_ARRAY: delete [] (real_t*)value; break;
                    case SETTING_TYPE_STRING: delete (std::string*)value; break;

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
                case SETTING_TYPE_STRING: return "a string";

                default: throw SettingsException("Unrecognized setting type: %d.", s);
            }
        }

        template<typename T>
        void _DefineSetting(const std::string& name, const std::string& desc, const T& defaultValue, enum setting_type type, bool mandatory);

        template<typename T>
        void _DefineSetting(const std::string& name, const std::string& desc, const len_t ndims, const len_t dims[], const T *defaultValue, enum setting_type type, bool mandatory);

        template<typename T>
        T *_GetArray(const std::string&, const len_t, len_t[], enum setting_type, bool markused=true);
        setting_t *_GetSetting(const std::string&, enum setting_type, bool markused=true);

        template<typename T>
        void _SetSetting(const std::string& name, const T& value, enum setting_type type);

        template<typename T>
        void _SetSetting(const std::string& name, const len_t ndims, const len_t dims[], T *value, enum setting_type type);

    public:
        Settings();
        ~Settings();

        const std::map<std::string, setting_t*> GetSettings() const { return this->settings; }

        // DEFINE SETTINGS
        void DefineSetting(const std::string& name, const std::string& desc, bool defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, int_t defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, len_t n, const int_t *defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, len_t ndims, const len_t dims[], const int_t *defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, real_t defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, len_t n, const real_t *defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, len_t ndims, const len_t dims[], const real_t *defaultValue, bool mandatory=false);
        void DefineSetting(const std::string& name, const std::string& desc, const std::string& defaultValue, bool mandatory=false);

		void UndefineSetting(const std::string& name);

        enum setting_type GetType(const std::string&);

        // GETTERS
        bool GetBool(const std::string&, bool markused=true);
        int_t GetInteger(const std::string&, bool markused=true);
        const int_t *GetIntegerArray(const std::string& name, const len_t nExpectedDims, len_t dims[], bool markused=true);
        real_t GetReal(const std::string&, bool markused=true);
        const real_t *GetRealArray(const std::string& name, const len_t nExpectedDims, len_t dims[], bool markused=true);
        const std::string GetString(const std::string&, bool markused=true);
        std::vector<std::string> GetStringList(const std::string&, const char delim=';', bool markused=true);

        void MarkUsed(const std::string&);

        // SETTERS
        void SetSetting(const std::string& name, bool value);
        void SetSetting(const std::string& name, int_t value);
        void SetSetting(const std::string& name, len_t n, int_t *value);
        void SetSetting(const std::string& name, len_t ndims, const len_t dims[], int_t *value);
        void SetSetting(const std::string& name, real_t value);
        void SetSetting(const std::string& name, len_t n, real_t *value);
        void SetSetting(const std::string& name, len_t ndims, const len_t dims[], real_t *value);
        void SetSetting(const std::string& name, const std::string& value);

        void DisplaySettings();
    };
}

#endif/*_DREAM_SETTINGS_HPP*/
