#ifndef _DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP
#define _DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class HotTailCurrentDensityFromDistributionFunction : public FVM::EquationTerm {
    private:
        FVM::Grid *fluidGrid;
        FVM::Grid *hottailGrid;
        
        FVM::UnknownQuantityHandler *unknowns;
        PitchScatterFrequency *nuD;
        real_t **nuD_vec = nullptr;

        len_t id_fhot;
        len_t id_Eterm;
        len_t id_ncold;
        len_t id_Tcold;
        len_t id_ni;
        
        real_t 
            **J1Weights = nullptr,
            **J2Weights,
            **diffWeights,
            **dNuDmat;
        
        real_t 
            *dEterm,
            *j1Vec,
            *j2Vec;

        len_t *np;

        bool hasBeenInitialised = false;
        bool isCollFreqModeFULL;

        void Deallocate();
        void SetJ1Weights(const real_t *Eterm, const real_t *const*nu_D, real_t **weights);
        void AddJacobianBlockMaxwellian(const len_t derivId, FVM::Matrix *jac);
    public:
        HotTailCurrentDensityFromDistributionFunction(
            FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, 
            FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD,
            enum OptionConstants::collqty_collfreq_mode collfreq_mode, bool
        );
        virtual ~HotTailCurrentDensityFromDistributionFunction();

        virtual len_t GetNumberOfNonZerosPerRow() const override
            { return hottailGrid->GetMomentumGrid(0)->GetNp1(); }

        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP*/
