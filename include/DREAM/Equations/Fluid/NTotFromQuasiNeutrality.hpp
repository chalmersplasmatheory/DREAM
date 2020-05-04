#ifndef _DREAM_N_TOT_FROM_QUASI_NEUTRALITY_HPP
#define _DREAM_N_TOT_FROM_QUASI_NEUTRALITY_HPP

#include "DREAM/IonHandler.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class NTotFromQuasiNeutrality : public FVM::PredeterminedParameter {
    private:
        IonHandler *ions;

    public:
        NTotFromQuasiNeutrality(FVM::Grid*, IonHandler*);
        virtual ~NTotFromQuasiNeutrality();

        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_N_TOT_FROM_QUASI_NEUTRALITY_HPP*/
