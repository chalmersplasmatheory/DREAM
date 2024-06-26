#ifndef _DREAM_FVM_INTERPOLATOR_1D_HPP
#define _DREAM_FVM_INTERPOLATOR_1D_HPP

#include "FVM/config.h"
#include "FVM/FVMException.hpp"

namespace DREAM::FVM {
    class Interpolator1D {
    public:
        enum interp_method {
            INTERP_NEAREST,
            INTERP_LINEAR,
            INTERP_LOGARITHMIC
        };

    private:
        len_t nx;       // Number of x points
        len_t nblocks;  // Number of blocks of nx points
                        // (i.e. # of y points is nx*nblocks)

        const real_t *x;
        const real_t *y;
        real_t *logy=nullptr;

        real_t *buffer=nullptr;

        bool xIncreasing;

        len_t _find_x(const real_t);
        const real_t *_eval_linear(const real_t);
        const real_t *_eval_logarithmic(const real_t);
        const real_t *_eval_nearest(const real_t);

        enum interp_method method;
		bool owns_data;

    public:
        Interpolator1D(const len_t, const len_t, const real_t*, const real_t*, enum interp_method meth=INTERP_LINEAR, bool owns_data=true);
        ~Interpolator1D();

        const real_t *Eval(const real_t);

        const real_t *GetBuffer() const { return this->buffer; }
        const real_t *GetX() const { return this->x; }
        const real_t *GetY() const { return this->y; }
    };

    class Interpolator1DException : public FVMException {
    public:
        template<typename ... Args>
        Interpolator1DException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("Interpolator1D");
        }
    };
}

#endif/*_DREAM_FVM_INTERPOLATOR_1D_HPP*/
