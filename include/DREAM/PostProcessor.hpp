#ifndef _DREAM_POST_PROCESSOR_HPP
#define _DREAM_POST_PROCESSOR_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"


namespace DREAM {
    class PostProcessor {
    private:
        FVM::Grid *fluidGrid;
        FVM::UnknownQuantityHandler *unknowns;

        // Arrays
        real_t *runawayRate;

        // IDs
        len_t id_n_re;

    public:
        PostProcessor(FVM::Grid*, FVM::UnknownQuantityHandler*);
        ~PostProcessor();

        // Getters
        const real_t *GetRunawayRate() const { return this->runawayRate; }
        
        void Process(const real_t t);
    };
}

#endif/*_DREAM_POST_PROCESSOR_HPP*/
