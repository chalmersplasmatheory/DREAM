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
        len_t id_jhot;
        len_t id_pcut;
        len_t id_Eterm;
        len_t id_ncold;
        len_t id_Tcold;
        len_t id_ni;
        

        len_t nr;
        len_t *np = nullptr;

        real_t **p;
        real_t **Delta_p;
        real_t **delta_p = nullptr;

        real_t **hWeights;
        real_t **gWeights;

        bool **useLorentzLimit;

        real_t **diffWeights;
        
        real_t *dEterm;
        real_t **dNuDmat;

        bool hasBeenInitialised = false;

        void Deallocate();
        void SetGWeights(const real_t *Eterm, const real_t *const*nu_D, real_t **weights);
    public:
        HotTailCurrentDensityFromDistributionFunction(
            FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, 
            FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD
        );
        virtual ~HotTailCurrentDensityFromDistributionFunction();

        virtual len_t GetNumberOfNonZerosPerRow() const 
            { return hottailGrid->GetMomentumGrid(0)->GetNp1(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const 
            { 
                return GetNumberOfNonZerosPerRow() /* fhot */ 
                + 1 /* Eterm */ + 1 /* ncold */ 
                + unknowns->GetUnknown(id_ni)->NumberOfMultiples() /* ni */ ; }

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_DENSITY_FROM_DISTRIBUTION_FUNCTION_HPP*/
