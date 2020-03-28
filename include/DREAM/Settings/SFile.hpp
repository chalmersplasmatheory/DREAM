#ifndef _DREAM_SFILE_HPP
#define _DREAM_SFILE_HPP

#include <string>
#include <softlib/SFile.h>
#include "FVM/config.h"
#include "DREAM/Settings/Settings.hpp"

namespace DREAM {
    class SettingsSFile {
    public:
        static void LoadSettings(Settings*, const std::string&);
        static void LoadSettings(Settings*, SFile*);

        static void LoadSetting(const std::string&, enum Settings::setting_type, const len_t, SFile*, Settings*);

        static void LoadBool(const std::string&, SFile*, Settings*);
        static void LoadInteger(const std::string&, SFile*, Settings*);
        static void LoadReal(const std::string&, SFile*, Settings*);
        static void LoadIntegerArray(const std::string&, const len_t, SFile*, Settings*);
        static void LoadRealArray(const std::string&, const len_t, SFile*, Settings*);
    };
}

#endif/*_DREAM_SFILE_HPP*/
