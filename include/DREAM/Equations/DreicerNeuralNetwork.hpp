#ifndef _DREAM_EQUATIONS_DREICER_NEURAL_NETWORK_HPP
#define _DREAM_EQUATIONS_DREICER_NEURAL_NETWORK_HPP

namespace DREAM { class DreicerNeuralNetwork; }

#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/config.h"

namespace DREAM {
    class DreicerNeuralNetwork {
    private:
        RunawayFluid *REFluid;

        // DATA
        static const real_t input_std[8], input_mean[8];
        static const real_t output_std[1], output_mean[1];
        static const real_t b1[20],    b2[20],     b3[20],     b4[20],     b5[1];
        static const real_t W1[20*8], W2[20*20], W3[20*20], W4[20*20], W5[1*20];

        void nn_layer(
            const len_t, const len_t,
            const real_t*, const real_t*,
            const real_t*, real_t*, bool applyTransferFunction=true
        );

    public:
        DreicerNeuralNetwork(RunawayFluid*);

        real_t RunawayRate(const len_t, const real_t, const real_t, const real_t);
        real_t RunawayRate_derived_params(
            const real_t, const real_t, const real_t,
            const real_t, const real_t, const real_t,
            const real_t, const real_t
        );
        real_t nn_fix(
            const real_t EDD, const real_t logTheta, 
            const real_t Zeff, const real_t Zeff0, const real_t Z0Z, 
            const real_t Z0_Z, const real_t logNfree, const real_t nfree_ntot,
            const real_t nfree, const real_t tauEE, bool use_fix); 

        bool IsApplicable(const real_t);
    };
}

#endif/*_DREAM_EQUATIONS_DREICER_NEURAL_NETWORK_HPP*/
