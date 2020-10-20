#ifndef _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP
#define _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/NIST.hpp"

namespace DREAM {
    class RadiatedPowerTerm : public FVM::DiagonalComplexTerm {
    private:
        ADAS *adas;
        NIST *nist;
        IonHandler *ionHandler;

        len_t id_ncold;
        len_t id_Tcold;
        len_t id_ni;
        
        bool with_PRB = true;
    protected:
        virtual len_t GetNumberOfWeightsElements() override 
            {return ionHandler->GetNzs() * grid->GetNCells();}

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        RadiatedPowerTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, ADAS*, NIST*);
    };
}


#endif /*_DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP*/
