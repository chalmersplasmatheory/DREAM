#ifndef _DREAM_IO_HPP
#define _DREAM_IO_HPP

#include <string>
#include <vector>
#include "DREAM/Simulation.hpp"

namespace DREAM {
    class IO {
    public:
        typedef enum {
            MESSAGE_GENERAL,

            WARNING_DREICER_NEURAL_NETWORK_INVALID,
            WARNING_INCOMPATIBLE_TRANSPORT,
            WARNING_KINETIC_AVALANCHE_NO_HOT_GRID,
            WARNING_NO_DEFAULT_PRECONDITIONING,
            WARNING_OVERRIDE_CUSTOM_P_GRID,
            WARNING_P_UNDERRESOLVED,
            WARNING_INCONSISTENT_RE_TRANSPORT,
            WARNING_TRITIUM_GENERATION_INVALID,
            WARNING_KINETIC_COMPTON_NO_HOT_GRID,
            WARNING_KINETIC_TRITIUM_NO_HOT_GRID,
			WARNING_RESIDUAL_NOT_CONVERGED,

            // [Must always come last]
            MESSAGE_LAST
        } message_t;

		static std::vector<std::string> emitted_warning_messages;
		static Simulation *simulation;

		static void Init(Simulation*);
		static void Deinit();

        static bool VerifyMessage(const message_t);

        static void PrintError(const char*);
        static void PrintError(const message_t, const char*);
        template<typename ... Args>
        static void PrintError(const message_t, const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 2, 0)))
        #endif
        ;
        template<typename ... Args>
        static void PrintError(const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 1, 0)))
        #endif
        ;

        static void PrintWarning(const char*);
        static void PrintWarning(const message_t, const char*);
        template<typename ... Args>
        static void PrintWarning(const message_t, const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 2, 0)))
        #endif
        ;
        template<typename ... Args>
        static void PrintWarning(const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 1, 0)))
        #endif
        ;

        static void PrintInfo();
        static void PrintInfo(const char*);
        static void PrintInfo(const message_t, const char*);
        template<typename ... Args>
        static void PrintInfo(const message_t, const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 2, 0)))
        #endif
        ;
        template<typename ... Args>
        static void PrintInfo(const char*, Args&& ...)
        #if defined(__GNUC__) || defined(__clang__)
        __attribute__((format (printf, 1, 0)))
        #endif
        ;

        static const std::string PRINT_YES, PRINT_NO;

    private:
        static bool message_checklist[MESSAGE_LAST];

    };
}

#include "IO.tcc"

#endif/*_DREAM_IO_HPP*/
