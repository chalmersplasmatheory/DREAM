#ifndef _DREAM_FVM_GSL_EXCEPTION_HPP
#define _DREAM_FVM_GSL_EXCEPTION_HPP

#include <string>
#include "FVM/FVMException.hpp"

namespace DREAM::FVM {
    class GSLException : public FVMException {
    private:
        std::string reason, file;
        int line, gsl_errno;

    public:
        template<typename ... Args>
        GSLException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("GSL");
        }
        GSLException(const char *rsn, const char *file, int line, int gsl_errno)
            : FVMException("%s:%d: %s (%s).", file, line, rsn, GSLException::GetErrorCodeName(gsl_errno).c_str()) {

            this->reason = rsn;
            this->file = file;
            this->line = line;
            this->gsl_errno = gsl_errno;
        }

        const std::string GetErrorCodeName() {
            return GSLException::GetErrorCodeName(this->gsl_errno);
        }
        static const std::string GetErrorCodeName(const int gsl_errno) {
            switch (gsl_errno) {
                case GSL_SUCCESS: return "GSL_SUCCESS";
                case GSL_FAILURE: return "GSL_FAILURE";
                case GSL_CONTINUE: return "GSL_CONTINUE";
                case GSL_EDOM: return "GSL_EDOM";
                case GSL_ERANGE: return "GSL_ERANGE";
                case GSL_EFAULT: return "GSL_EFAULT";
                case GSL_EINVAL: return "GSL_EINVAL";
                case GSL_EFAILED: return "GSL_EFAILED";
                case GSL_EFACTOR: return "GSL_EFACTOR";
                case GSL_ESANITY: return "GSL_ESANITY";
                case GSL_ENOMEM: return "GSL_ENOMEM";
                case GSL_EBADFUNC: return "GSL_EBADFUNC";
                case GSL_ERUNAWAY: return "GSL_ERUNAWAY";
                case GSL_EMAXITER: return "GSL_EMAXITER";
                case GSL_EZERODIV: return "GSL_EZERODIV";
                case GSL_EBADTOL: return "GSL_EBADTOL";
                case GSL_ETOL: return "GSL_ETOL";
                case GSL_EUNDRFLW: return "GSL_EUNDRFLW";
                case GSL_EOVRFLW: return "GSL_EOVRFLW";
                case GSL_ELOSS: return "GSL_ELOSS";
                case GSL_EROUND: return "GSL_EROUND";
                case GSL_EBADLEN: return "GSL_EBADLEN";
                case GSL_ENOTSQR: return "GSL_ENOTSQR";
                case GSL_ESING: return "GSL_ESING";
                case GSL_EDIVERGE: return "GSL_EDIVERGE";
                case GSL_EUNSUP: return "GSL_EUNSUP";
                case GSL_EUNIMPL: return "GSL_EUNIMPL";
                case GSL_ECACHE: return "GSL_ECACHE";
                case GSL_ETABLE: return "GSL_ETABLE";
                case GSL_ENOPROG: return "GSL_ENOPROG";
                case GSL_ENOPROGJ: return "GSL_ENOPROGJ";
                case GSL_ETOLF: return "GSL_ETOLF";
                case GSL_ETOLX: return "GSL_ETOLX";
                case GSL_ETOLG: return "GSL_ETOLG";
                case GSL_EOF: return "GSL_EOF";
                
                default:
                    return "UNRECOGNIZED CODE";
            }
        }
    };
}

#endif/*_DREAM_FVM_GSL_EXCEPTION_HPP*/
