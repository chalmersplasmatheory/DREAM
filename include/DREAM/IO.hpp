#ifndef _DREAM_IO_HPP
#define _DREAM_IO_HPP

#include <string>

namespace DREAM {
    class IO {
    public:
        typedef enum {
            MESSAGE_GENERAL,
            // [Must always come last]
            MESSAGE_LAST
        } message_t;

        static bool VerifyMessage(const message_t);

        template<typename ... Args>
        static void PrintError(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintError(const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintWarning(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintWarning(const std::string&, Args&& ...);
        static void PrintInfo();
        template<typename ... Args>
        static void PrintInfo(const message_t, const std::string&, Args&& ...);
        template<typename ... Args>
        static void PrintInfo(const std::string&, Args&& ...);

        static const std::string PRINT_YES, PRINT_NO;

    private:
        static bool message_checklist[MESSAGE_LAST];

    };
}

#include "IO.tcc"

#endif/*_DREAM_IO_HPP*/
