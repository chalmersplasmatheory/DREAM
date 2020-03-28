/**
 * Implementation of some screen I/O helper routines.
 */

#include <omp.h>
#include "DREAM/IO.hpp"


using namespace DREAM;

bool IO::message_checklist[MESSAGE_LAST] = {false};
#ifdef COLOR_TERMINAL
    const std::string IO::PRINT_YES = "\x1B[1;32mYES\x1B[0m";
    const std::string IO::PRINT_NO  = "\x1B[1;32mNO\x1B[0m";
#else
    const std::string IO::PRINT_YES = "YES";
    const std::string IO::PRINT_NO  = "NO";
#endif

/**
 * Print a single new-line character in the 'Info' channel.
 */
void IO::PrintInfo() {
    PrintInfo(MESSAGE_GENERAL, "");
}

/**
 * Verify the given message ID against the checklist
 * that the message has not previously been emitted
 * during the run.
 *
 * id: Message ID.
 *
 * RETURNS true if the message has NOT been emitted
 * previously (and is safe to emit now).
 */
bool IO::VerifyMessage(const message_t id) {
    bool retval = false,
         valid = (MESSAGE_GENERAL < id && id < MESSAGE_LAST);

    if (!valid)
        return true;
    else if (message_checklist[id])
        return false;

    #pragma omp critical (DREAM_VerifyMessage)
    {
        if (message_checklist[id] == false) {
            message_checklist[id] = true;
            retval = true;
        }
    }

    return retval;
}

