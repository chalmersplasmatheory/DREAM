#ifndef _DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_PCUT_TERM_HPP
#define _DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_PCUT_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Equations/PitchScatterFrequency.hpp"
#include "DREAM/Constants.hpp"

namespace DREAM {
    class HotTailCurrentPCutTerm : public FVM::EquationTerm {
    private:
        FVM::Grid *fluidGrid;
        FVM::Grid *hottailGrid;
        
        FVM::UnknownQuantityHandler *unknowns;
        PitchScatterFrequency *nuD;

        len_t id_fhot;
        len_t id_jhot;
        len_t id_pcut;
        len_t id_Eterm;
        len_t id_ncold;
        len_t id_ni;

        len_t nr;
        real_t *dp = nullptr;

        real_t *GOverH;
//        real_t *dGOverH;
        len_t *iCut;

        bool hasBeenInitialised = false;

        real_t EvalGOverH(len_t ir, len_t derivId, len_t n=0);
        void Deallocate();

    public:
        HotTailCurrentPCutTerm(
            FVM::Grid *fluidGrid, FVM::Grid *hottailGrid, 
            FVM::UnknownQuantityHandler *u, PitchScatterFrequency *nuD
        );
        virtual ~HotTailCurrentPCutTerm();

        virtual len_t GetNumberOfNonZerosPerRow() const { return 3; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const 
        { 
            return GetNumberOfNonZerosPerRow() /* fhot */ 
            + 1 /* Eterm */ + 1 /* ncold */ + 1 /* pcut */ 
            + unknowns->GetUnknown(id_ni)->NumberOfMultiples() /* ni */ ; 
        }

        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

    };
}

#endif/*_DREAM_EQUATION_FLUID_HOT_TAIL_CURRENT_PCUT_TERM_HPP*/
