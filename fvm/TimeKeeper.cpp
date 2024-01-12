/**
 * The TimeKeeper class keeps a number of DurationTimer's which can be started
 * and stopped during a simulation. In the end, the timing information can be
 * printed in a formatted way, and stored to an SFile object.
 */

#include <softlib/SFile.h>
#include "FVM/TimeKeeper.hpp"


using namespace DREAM::FVM;
using namespace std;


/**
 * Constructor.
 *
 * name: Name of TimeKeeper (only used when printing information to stdout).
 */
TimeKeeper::TimeKeeper(const string& name)
    : name(name) {}

/**
 * Destructor.
 */
TimeKeeper::~TimeKeeper() {
    for (auto it : timers)
        delete it;
}

/**
 * Add a new timer.
 *
 * shortname: Short name of new timer (used for saving to file).
 * longname:  A slightly longer name of the new timer (used for printing to stdout).
 *
 * RETURNS the ID of the new timer.
 */
len_t TimeKeeper::AddTimer(const string& shortname, const string& longname) {
    timers.push_back(new (struct tk){new DurationTimer(), longname, shortname});
    return timers.size()-1;
}

/**
 * Reset the specified timer.
 *
 * timer: ID of timer to reset.
 */
void TimeKeeper::ResetTimer(const len_t timer) {
    timers[timer]->timer->Reset();
}

/**
 * Start the specified timer.
 *
 * timer: ID of timer to start.
 */
void TimeKeeper::StartTimer(const len_t timer) {
    timers[timer]->timer->Start();
}

/**
 * Stop the specified timer.
 *
 * timer: ID of the timer to stop.
 */
void TimeKeeper::StopTimer(const len_t timer) {
    timers[timer]->timer->Stop();
}

/**
 * Print the results accumulated by this timer.
 *
 * printTitle:  If true, prints a title for the information.
 * normalizeto: If non-negative, prints all timings as normalized to the
 *              timer with the specified ID.
 */
void TimeKeeper::PrintTimings(bool printTitle, const int_t normalizeto) {
    len_t timer=0;
    bool normalize = (normalizeto >= 0);

    if (normalize)
        timer = static_cast<len_t>(normalizeto);

    // Determine longest name of timers
    int maxlen = 0;
    for (auto s : timers)
        maxlen = max(static_cast<int>(s->longname.length()), maxlen);

    real_t nrm = (normalize ? timers[timer]->timer->GetMicroseconds() : 1);

    if (printTitle)
        printf("[Timing for %s]\n", this->name.c_str());

    real_t sum = 0;
    for (auto tm : timers) {
        DurationTimer *dt = tm->timer;
        real_t s = dt->GetMicroseconds();
        
        // Print timing
        if (normalize) {
            if (s == nrm) continue;

            sum += s;
            printf("  %-*s  %3.2f%%\n", maxlen+1, (tm->longname+":").c_str(), s/nrm*100);
        } else {
            real_t ms = s/1000;

            if (ms > 1000)
                printf("  %-*s  %3.4f s", maxlen+1, (tm->longname+":").c_str(), ms/1000);
            else
                printf("  %-*s  %3.4f ms", maxlen+1, (tm->longname+":").c_str(), ms);
        }
    }

    // Print residual (in normalized mode; otherwise people should easily
    // be able to figure out the residual themselves; arithmetic is healthy)
    if (normalize)
        printf("  %-*s  %3.2f%%\n", maxlen+1, "Other work:", (nrm-sum)/nrm*100);
}

/**
 * Save timing information to the given SFile object.
 *
 * sf:   SFile object to save timing information to.
 * path: Path in SFile object to save information to.
 */
void TimeKeeper::SaveTimings(SFile *sf, const std::string& path) {
    for (auto tm : timers) {
        DurationTimer *dt = tm->timer;
        std::string dsetname = path+"/"+tm->shortname;

        sf->WriteScalar(dsetname, dt->GetMicroseconds());
        sf->WriteAttribute_string(dsetname, "desc", tm->longname);
    }
}

