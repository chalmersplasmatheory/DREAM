#ifndef _DREAM_EXCEPTION_HPP
#define _DREAM_EXCEPTION_HPP

#include "FVM/FVMException.hpp"

namespace DREAM {
    class DREAMException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        DREAMException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("DREAM exception");
        }
    };
}

#endif/*_DREAM_EXCEPTION_HPP*/
