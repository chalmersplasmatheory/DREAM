#ifndef _DREAM_FVM_TIME_KEEPER_HPP
#define _DREAM_FVM_TIME_KEEPER_HPP

#include <softlib/SFile.h>
#include <string>
#include <vector>
#include "FVM/DurationTimer.hpp"

namespace DREAM::FVM {
    class TimeKeeper {
    protected:
        struct tk {
            DurationTimer *timer;
            std::string longname;
            std::string shortname;

            ~tk() { delete timer; }
        };
        std::string name;
        std::vector<struct tk*> timers;

    public:
        TimeKeeper(const std::string&);
        ~TimeKeeper();

        len_t AddTimer(const std::string&, const std::string&);
        void ResetTimer(const len_t);
        void StartTimer(const len_t);
        void StopTimer(const len_t);

        void PrintTimings(bool printTitle=true, const int_t normalizeto=0);
        void SaveTimings(SFile*, const std::string& path="");
    };
}

#endif/*_DREAM_FVM_TIME_KEEPER_HPP*/
