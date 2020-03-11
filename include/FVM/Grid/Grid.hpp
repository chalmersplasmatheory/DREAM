#ifndef _TQS_FVM_GRID_HPP
#define _TQS_FVM_GRID_HPP

#include <string>
#include "config.h"


namespace TQS::FVM {
    template<len_t N>
    class Grid {
    protected:
        // Convert index tuple to linear index
        template<typename ... Args>
        virtual len_t get_index(Args&& ..., len_t) const = 0;

        // Private, abstract methods
        virtual real_t get_xn(const len_t, Args&& ...) const = 0;
        virtual real_t get_xn_f(const len_t, Args&& ...) const = 0;

    public:
        // Copy constructor
        Grid(const len_t[N]);
        virtual ~Grid() {}

        template<typename ... Args>
        real_t GetX(const len_t, Args&& ...) const;
        template<typename ... Args>
        real_t GetX_f(const len_t, Args&& ...) const;

        virtual bool Rebuild(const real_t);

        template<typename ... Args>
        virtual real_t TQS::FVM::Volume(Args&& ...) const = 0;
    };
}

#include "FVM/Grid/Grid.tcc"

#endif/*_TQF_FVM_GRID_HPP*/
