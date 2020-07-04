#ifndef _DREAM_FVM_DURATION_TIMER_HPP
#define _DREAM_FVM_DURATION_TIMER_HPP

#include <chrono>
#include <string>
#include "FVM/config.h"

namespace DREAM::FVM {
    class DurationTimer {
    protected:
        std::chrono::time_point<std::chrono::high_resolution_clock> tic;
        std::chrono::microseconds accum;

    public:
        DurationTimer();

        real_t GetMicroseconds() const;
        real_t GetMilliseconds() const;

        void Reset();
        void Start(bool reset=false);
        void Stop();
    };
}
#endif/*_DREAM_FVM_DURATION_TIMER_HPP*/
