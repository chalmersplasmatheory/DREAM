#ifndef _DREAM_SETTINGS_LOAD_DATA_HPP
#define _DREAM_SETTINGS_LOAD_DATA_HPP

#include "FVM/Interpolator1D.hpp"
#include "FVM/Interpolator3D.hpp"

namespace DREAM {
    struct dream_2d_data {
        real_t
            *x, *t, *r;
        len_t nt, nr;
        enum FVM::Interpolator1D::interp_method interp;
    };
    struct dream_4d_data {
        real_t
            **x,        // Size nt-by-(nr*np1*np2)
            *t, *r, *p1, *p2;
        len_t nt, nr, np1, np2;
        enum FVM::Interpolator3D::momentumgrid_type gridtype;
        enum FVM::Interpolator3D::interp_method ps_interp;
        enum FVM::Interpolator1D::interp_method time_interp;
    };
}

#endif/*_DREAM_SETTINGS_LOAD_DATA_HPP*/
