#ifndef _DREAM_ADAS_HPP
#define _DREAM_ADAS_HPP

#include <unordered_map>
#include <gsl/gsl_interp.h>
#include "DREAM/ADASRateInterpolator.hpp"
#include "FVM/FVMException.hpp"

namespace DREAM {
    class ADAS {
    private:
        std::unordered_map<len_t, ADASRateInterpolator**> intp;

        static const len_t
            IDX_ACD, IDX_CCD, IDX_SCD, IDX_PLT, IDX_PRB;

        len_t get_isotope_index(const len_t, const len_t A=0) const;
        std::unordered_map<len_t, ADASRateInterpolator**>::const_iterator
            get_element(const len_t, const len_t A=0) const;

        // Maximum possible atomic mass (or larger; should just be an
        // arbitrary large number to help with indexing...)
        static const len_t MAX_ATOMIC_MASS = 250;

    public:
        ADAS(const gsl_interp2d_type *interp=gsl_interp2d_bicubic);
        ~ADAS();

        bool HasElement(const len_t Z, const len_t A=0) const;

        ADASRateInterpolator *GetACD(const len_t Z, const len_t A=0) const;
        ADASRateInterpolator *GetCCD(const len_t Z, const len_t A=0) const;
        ADASRateInterpolator *GetSCD(const len_t Z, const len_t A=0) const;
        ADASRateInterpolator *GetPLT(const len_t Z, const len_t A=0) const;
        ADASRateInterpolator *GetPRB(const len_t Z, const len_t A=0) const;

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
