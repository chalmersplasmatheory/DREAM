#ifndef _DREAM_EQUATION_FLUID_EXTERNAL_AVALANCHE_TERM_HPP
#define _DREAM_EQUATION_FLUID_EXTERNAL_AVALANCHE_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceRP.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
namespace DREAM {
    class ExternalAvalancheTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        const real_t pCutoff, aExp, scaleFactor;
        len_t id_ntot;
        len_t id_Efield;

        real_t *dGammaFluid = nullptr;
        len_t nr_tmp=0;
        void AllocateDGamma(){ // if nr has changed, (re)allocate dGamma
            if(nr_tmp != this->nr){
                if(dGammaFluid != nullptr)
                    delete [] dGammaFluid;
                dGammaFluid = new real_t[this->nr];
                nr_tmp = this->nr;
            }
        }
    protected:
        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;
    public:
        ExternalAvalancheTerm(FVM::Grid*, real_t pc, real_t aExp, RunawayFluid*, FVM::UnknownQuantityHandler*, real_t scaleFactor=1.0);
		~ExternalAvalancheTerm();
    };
}


#endif /*_DREAM_EQUATION_FLUID_EXTERNAL_AVALANCHE_TERM_HPP*/
