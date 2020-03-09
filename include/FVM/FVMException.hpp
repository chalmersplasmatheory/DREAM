#ifndef _TQS_FVM_EXCEPTION_H
#define _TQS_FVM_EXCEPTION_H

#include <string>
#include <vector>
#include <softlib/SOFTLibException.h>

/**
 * General-purpose TQS::FVM Exception class.
 */
namespace TQS::FVM {
    class FVMException : public SOFTLibException {
    private:
        std::vector<std::string> modules;
    public:
        FVMException() {}

        template<typename ... Args>
        FVMException(const std::string& msg, Args&& ... args)
            : SOFTLibException(msg, std::forward<Args>(args)...) {}

        void AddModule(const std::string& m) { modules.push_back(m); }

        std::vector<std::string> &GetModules() { return modules; }
        std::string GetModulesString() {
            std::string m;
            if (modules.size() > 0) {
                std::string m = modules.front();
                for (std::vector<std::string>::iterator it = modules.begin()+1; it != modules.end(); it++)
                    m += "/" + (*it);
            }

            return m;
        }
    };
}

#endif/*_TQS_FVM_EXCEPTION_H*/
