#ifndef _DREAM_EQUATION_FLUID_LCFS_LOSS_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_LCFS_LOSS_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

namespace DREAM {
    class LCFSLossRateTerm : public RunawaySourceTerm, public FVM::DiagonalComplexTerm {
    private:
    
        FVM::UnknownQuantityHandler *unknowns;
        real_t scaleFactor;
        real_t* t_loss;
        const len_t userGivenPsiEdge_t0;
        real_t psi_edge_t0;
        const len_t id_psi;
        FVM::RadialGrid *rGrid;
        real_t* GammaLoss = nullptr;
        int_t ir_LCFS;
        bool signFixed = false;
        real_t sign = -1;
        
        real_t PsiDiff(len_t);
        real_t InterpolatePsi(len_t);
        void FindRadiusOfLCFS();
        real_t StepFunction(len_t);
        void SetGammaLoss();
        void Deallocate();
        
        /***************************** FUNCTIONS NEEDED FOR DIAGONALCOMPLEXTERM *****************************/
        virtual void SetWeights() override; 
        
        // Now empty but could be made to include the psi dependence:
        virtual void SetDiffWeights(len_t /*derivID*/, len_t /*nMultiples*/) override; 
        /****************************************************************************************************/
        
        
    public:
    
        LCFSLossRateTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, FVM::Grid* /*operandGrid*/, real_t, real_t*, len_t userGivenPsiEdge_t0=0, real_t PsiEdge_t0=0);
        ~LCFSLossRateTerm();
        // USES REBUILD, SETVECTORELEMENT, SETJACOBIANBLOCK FROM PARENT CLASSES
        
        virtual bool GridRebuilt() override;
        virtual const real_t *GetLCFSLossWeights();
    
    };
}

#endif/*_DREAM_EQUATION_FLUID_LCFS_LOSS_RATE_TERM_HPP*/
