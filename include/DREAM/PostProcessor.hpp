#ifndef _DREAM_POST_PROCESSOR_HPP
#define _DREAM_POST_PROCESSOR_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "FVM/Equation/MomentQuantity.hpp"

namespace DREAM {
    class PostProcessor {
    private:
        FVM::Grid *fluidGrid;
        FVM::UnknownQuantityHandler *unknowns;

        // Arrays
        real_t *runawayRate;

        // hot-electron region definition
        real_t pThreshold;
        FVM::MomentQuantity::pThresholdMode pThresholdMode;

        // IDs
        len_t id_n_re;

    public:
        PostProcessor(
            FVM::Grid*, FVM::UnknownQuantityHandler*, real_t pThreshold, 
            FVM::MomentQuantity::pThresholdMode
        );
        ~PostProcessor();

        // Getters
        const real_t *GetRunawayRate() const { return this->runawayRate; }
        const real_t GetPThreshold() const {return pThreshold;}
        const FVM::MomentQuantity::pThresholdMode GetPThresholdMode() const {return pThresholdMode;}
        
        void Process(const real_t t);
    };
}

#endif/*_DREAM_POST_PROCESSOR_HPP*/
