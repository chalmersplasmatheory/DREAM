
#include <string>
#include <vector>
#include "FVM/config.h"

namespace DREAM {

/**
 * Prints an error message to stderr.
 * Supports printf syntax.
 *
 * id:  Optional error ID (to only print error once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void IO::PrintError(const char *msg, Args&& ... args) {
    PrintError(IO::MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}

template<typename ... Args>
void IO::PrintError(const IO::message_t id, const char *msg, Args&& ... args) {
    if (!IO::VerifyMessage(id))
        return;

#ifdef COLOR_TERMINAL
    fputs("\x1B[1;31m[ERROR]\x1B[0m ", stderr);
#else
    fputs("[ERROR] ", stderr);
#endif
    fprintf(stderr, msg, std::forward<Args>(args) ...);
    fputc('\n', stderr);
}

/**
 * Prints a warning message once (if the same warning
 * has already been printed, this function doesn't
 * print anything).
 * Supports printf syntax.
 *
 * id:  Optional warning ID (to only print warning once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void IO::PrintWarning(const char *msg, Args&& ... args) {
    PrintWarning(IO::MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}
template<typename ... Args>
void IO::PrintWarning(const IO::message_t id, const char *msg, Args&& ... args) {
    if (!IO::VerifyMessage(id))
        return;

	// Generate warning string
	int size = std::snprintf(nullptr, 0, msg, args ...) + 1;
	if (size <= 0)
		return;
	char *buf = new char[size];
	std::snprintf(buf, size, msg, args ...);
	std::string smsg(buf);

	bool prepended = false;
	if (IO::simulation != nullptr) {
		EquationSystem *eqsys = IO::simulation->GetEquationSystem();
		if (eqsys != nullptr) {
			FVM::UnknownQuantityHandler *u = eqsys->GetUnknownHandler();
			if (u != nullptr && u->Size() > 0) {
				smsg = "At time step " +
					std::to_string(u->GetUnknown(0)->GetNumberOfSavedSteps()+1) +
					": " + smsg;
				prepended = true;
			}
		}
	}
	
	if (!prepended)
		smsg = "During initialization: " + smsg;

	emitted_warning_messages.push_back(smsg);

#ifdef COLOR_TERMINAL
    fputs("\x1B[1;33m[WARNING]\x1B[0m ", stderr);
#else
    fputs("[WARNING] ", stderr);
#endif
    fputs(smsg.c_str(), stderr);
    fputc('\n', stderr);
}

/**
 * Prints a general information message to stdout.
 * Supports printf syntax.
 *
 * id:  Optional message ID (to only print message once).
 * msg: Message (in printf-format) to print.
 */
template<typename ... Args>
void IO::PrintInfo(const char *msg, Args&& ... args) {
    PrintInfo(MESSAGE_GENERAL, msg, std::forward<Args>(args) ...);
}
template<typename ... Args>
void IO::PrintInfo(const IO::message_t id, const char *msg, Args&& ... args) {
    if (!IO::VerifyMessage(id))
        return;
    
    fprintf(stdout, msg, std::forward<Args>(args) ...);
    fputc('\n', stdout);
}

}//namespace DREAM
