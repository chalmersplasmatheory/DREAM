#ifndef _DREAM_EQUATION_FLUID_IONISATION_HEATING_TERM_HPP
#define _DREAM_EQUATION_FLUID_IONISATION_HEATING_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/NIST.hpp"

namespace DREAM {
    class IonisationHeatingTerm : public FVM::DiagonalComplexTerm {
    private:
        ADAS *adas;
        NIST *nist;
        IonHandler *ionHandler;

        bool initParamsInitialised = false;
        real_t *T_cold_init = nullptr;
        //real_t *n_cold_init = nullptr;
        len_t nr_init=0;

        len_t id_ncold;
        len_t id_nhot;
        len_t id_Tcold;
        len_t id_ni;

        void AllocateInitParams();
    protected:
        virtual len_t GetNumberOfWeightsElements() override 
            {return ionHandler->GetNzs() * grid->GetNCells();}

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        IonisationHeatingTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, ADAS*, NIST*);
    };
}


#endif /*_DREAM_EQUATION_FLUID_IONISATION_HEATING_TERM_HPP*/
