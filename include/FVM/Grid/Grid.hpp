#ifndef _TQS_FVM_GRID_HPP
#define _TQS_FVM_GRID_HPP

#include <string>
#include "config.h"
#include "FVM/Grid/Grid1D.hpp"


namespace TQS::FVM {
    template<int N>
    class Grid {
    private:
        // List of pointers to Grid1D objects
        Grid1D *dimensions[N];

        real_t *volumes;

        len_t get_index(len_t) const;
        template<typename ... Args>
        len_t get_index(Args&& ... args, len_t) const;

        void insert_dimensions(Grid1D*) const;
        template<typename ... Args>
        void insert_dimensions(Grid1D*, Args&& ... args) const;

    public:
        Grid(Grid*);

        template<typename ... Args>
        Grid(Args&& ... args) {
            this->insert_dimensions(args...);
        }

        bool Rebuild(const real_t);

        template<typename ... Args>
        real_t TQS::FVM::Volume(Args&& ... args) const;
    };
}

#endif/*_TQF_FVM_GRID_HPP*/
