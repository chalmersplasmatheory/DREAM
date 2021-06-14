#ifndef _DREAM_ADAS_HPP
#define _DREAM_ADAS_HPP

#include <map>
#include <gsl/gsl_interp.h>
#include "DREAM/ADASRateInterpolator.hpp"
#include "FVM/FVMException.hpp"

namespace DREAM {
    class ADAS {
    private:
        std::map<len_t, ADASRateInterpolator**> intp;

        static const len_t
            IDX_ACD, IDX_CCD, IDX_SCD, IDX_PLT, IDX_PRB;

        std::map<len_t, ADASRateInterpolator**>::const_iterator
            get_element(const len_t) const;

    public:
        ADAS(const gsl_interp2d_type *interp=gsl_interp2d_bicubic);
        ~ADAS();

        bool HasElement(const len_t Z) const;

        ADASRateInterpolator *GetACD(const len_t Z) const;
        ADASRateInterpolator *GetCCD(const len_t Z) const;
        ADASRateInterpolator *GetSCD(const len_t Z) const;
        ADASRateInterpolator *GetPLT(const len_t Z) const;
        ADASRateInterpolator *GetPRB(const len_t Z) const;

        void PrintElements() const;
    };

    class ADASException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        ADASException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("ADAS");
        }
    };
}

#endif/*_DREAM_ADAS_HPP*/
