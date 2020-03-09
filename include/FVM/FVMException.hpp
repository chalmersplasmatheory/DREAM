#ifndef _TQS_FVM_EXCEPTION_H
#define _TQS_FVM_EXCEPTION_H

#include <string>
#include <vector>
#include <exception>

/**
 * General-purpose TQS::FVM Exception class.
 */
namespace TQS::FVM {
    class FVMException : public std::exception {
    private:
        std::vector<std::string> modules;
    public:
        FVMException() {}
        FVMException(const std::string&);

		template<typename ... Args>
		FVMException(const std::string& msg, Args&& ... args) {
            ConstructErrorMessage(msg, std::forward<Args>(args) ...);
		}
        template<typename ... Args>
        void ConstructErrorMessage(const std::string& msg, Args&& ... args) {
			#define SOFTLIBEXCEPTION_BUFFER_SIZE 1000
			int n;
			char *buffer;

			buffer = (char*)malloc(sizeof(char)*SOFTLIBEXCEPTION_BUFFER_SIZE);

			n = snprintf(buffer, SOFTLIBEXCEPTION_BUFFER_SIZE, msg.c_str(), std::forward<Args>(args) ...);
			if (n < 0) {	// Encoding error. Just save error message unformatted.
				this->errormsg = msg;
			} else if (n >= SOFTLIBEXCEPTION_BUFFER_SIZE) {
				buffer = (char*)realloc(buffer, sizeof(char)*(n+1));
				n = snprintf(buffer, n, msg.c_str(), std::forward<Args>(args) ...);
				this->errormsg = buffer;
			} else this->errormsg = buffer;
        }

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

		virtual const char* what() const throw();
		virtual const std::string& whats() const throw();
	private:
		std::string errormsg;
    };
}

#endif/*_TQS_FVM_EXCEPTION_H*/
