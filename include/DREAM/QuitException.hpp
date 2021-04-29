#ifndef _DREAM_QUIT_EXCEPTION_HPP
#define _DREAM_QUIT_EXCEPTION_HPP

#include "FVM/FVMException.hpp"

namespace DREAM {
    class QuitException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        QuitException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("SIGQUIT");
        }
    };
}

#endif/*_DREAM_QUIT_EXCEPTION_HPP*/
