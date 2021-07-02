#ifndef _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP
#define _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/AMJUEL.hpp"

namespace DREAM {
    class RadiatedPowerTerm : public FVM::DiagonalComplexTerm {
    private:
        ADAS *adas;
        NIST *nist;
        AMJUEL *amjuel;
        IonHandler *ionHandler;
        
        enum OptionConstants::ion_opacity_mode *opacity_modes;

        len_t id_ncold;
        len_t id_Tcold;
        len_t id_ni;
        
        real_t 
            bremsPrefactor,
            bremsRel1,
            bremsRel2;
            
        bool includePRB = true;
    protected:
        virtual len_t GetNumberOfWeightsElements() override 
            {return ionHandler->GetNzs() * grid->GetNCells();}

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

        void SetWeights(const real_t*, real_t *w=nullptr);
        void SetDiffWeights(len_t derivId, len_t nMultiples, const real_t*);

    public:
        RadiatedPowerTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, ADAS*, NIST*, AMJUEL*,enum OptionConstants::ion_opacity_mode*, bool);
    };
}


#endif /*_DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP*/
