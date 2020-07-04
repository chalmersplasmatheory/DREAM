/**
 * Implementation of a timer that can be started and stopped several times,
 * gradually accumulating time.
 */

#include <chrono>
#include <string>
#include <sstream>
#include "FVM/DurationTimer.hpp"

using namespace std;
using namespace std::chrono;
using namespace std::chrono_literals;
using namespace DREAM::FVM;


/**
 * Constructor.
 */
DurationTimer::DurationTimer() {
    this->Reset();
}

/**
 * Get the number of microseconds which have been measured.
 */
real_t DurationTimer::GetMicroseconds() const {
    return static_cast<real_t>(this->accum.count());
}

/**
 * Get the number of milliseconds which have been measured.
 */
real_t DurationTimer::GetMilliseconds() const {
    return this->GetMicroseconds()*1e-3;
}

/**
 * Reset the timer.
 */
void DurationTimer::Reset() {
    this->accum = 0us;
}

/**
 * Start the timer.
 */
void DurationTimer::Start(bool reset) {
    if (reset)
        this->Reset();

    this->tic = high_resolution_clock::now();
}

/**
 * Stop the timer.
 */
void DurationTimer::Stop() {
    auto toc = high_resolution_clock::now();

    this->accum += duration_cast<microseconds>(toc - this->tic);
}

