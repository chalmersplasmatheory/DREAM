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
        static void LoadString(const std::string&, SFile*, Settings*);
        static void LoadIntegerArray(const std::string&, const len_t, SFile*, Settings*);
        static void LoadRealArray(const std::string&, const len_t, SFile*, Settings*);

        static void CreateGroup(const std::string&, const std::string&, std::vector<std::string>&, SFile*);
        static void SaveSettings(Settings*, SFile*, const std::string&);
        static void SaveBool(const std::string&, bool, SFile*);
        static void SaveInteger(const std::string&, int_t, SFile*);
        static void SaveReal(const std::string&, real_t, SFile*);
        static void SaveIntegerArray(const std::string&, const int_t*, const len_t, const len_t*, SFile*);
        static void SaveRealArray(const std::string&, const real_t*, const len_t, const len_t*, SFile*);
        static void SaveString(const std::string&, const std::string&, SFile*);
    };
}

#endif/*_DREAM_SFILE_HPP*/
