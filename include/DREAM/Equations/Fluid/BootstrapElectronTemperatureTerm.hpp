#ifndef _DREAM_EQUATION_FLUID_BOOTSTRAP_ELECTRON_TERMPERATURE_TERM_HPP
#define _DREAM_EQUATION_FLUID_BOOTSTRAP_ELECTRON_TERMPERATURE_TERM_HPP

#include "FVM/Grid/Grid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/Equations/BootstrapCurrent.hpp"
#include "DREAM/Equations/Fluid/BootstrapEquationTerm.hpp"

namespace DREAM {
    class BootstrapElectronTemperatureTerm : public BootstrapEquationTerm {

    private:
        BootstrapCurrent *Bootstrap;

    public:
        BootstrapElectronTemperatureTerm(FVM::Grid*, BootstrapCurrent*, IonHandler*, real_t sf=1.);

        virtual real_t GetCoefficient(len_t, len_t) override;
        virtual real_t GetPartialCoefficient(len_t, len_t, len_t, len_t) override;
    };
}

#endif /*_DREAM_EQUATION_FLUID_BOOTSTRAP_ELECTRON_TERMPERATURE_TERM_HPP*/
