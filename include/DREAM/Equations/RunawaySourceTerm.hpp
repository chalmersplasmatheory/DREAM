#ifndef _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP
#define _DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class RunawaySourceTerm {
    private:
        FVM::Grid *rst_grid;
        FVM::UnknownQuantityHandler *rst_unknowns;
        len_t id_E_field;

    public:
        RunawaySourceTerm(FVM::Grid*, FVM::UnknownQuantityHandler*);

        len_t GetXiIndexForEDirection(const len_t);
        real_t GetVolumeScaleFactor(const len_t);
    };
}

#endif/*_DREAM_EQUATIONS_RUNAWAY_SOURCE_TERM_HPP*/
