#ifndef _DREAM_NOT_IMPLEMENTED_EXCEPTION_HPP
#define _DREAM_NOT_IMPLEMENTED_EXCEPTION_HPP

#include "FVM/FVMException.hpp"

namespace DREAM {
    class NotImplementedException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        NotImplementedException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Not-yet-implemented");
        }
    };
}

#endif/*_DREAM_NOT_IMPLEMENTED_EXCEPTION_HPP*/
