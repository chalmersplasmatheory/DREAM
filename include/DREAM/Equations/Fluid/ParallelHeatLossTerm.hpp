#ifndef _DREAM_EQUATION_FLUID_HALO_REGION_HEAT_LOSS_TERM_HPP
#define _DREAM_EQUATION_FLUID_HALO_REGION_HEAT_LOSS_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTerm.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"

namespace DREAM {
    class HaloRegionHeatLossTerm : public FVM::DiagonalComplexTerm {
    private:
    
        FVM::UnknownQuantityHandler *unknowns;
        IonHandler *ions;

        len_t id_T_cold, id_W_cold, id_N_i, id_W_i, id_jtot;

        real_t scaleFactor;
        bool userGivenPsiEdge_t0;
        real_t psi_edge_t0;
        const len_t id_psi;
        FVM::RadialGrid *rGrid;
        real_t* HeatLoss = nullptr;
        int_t ir_LCFS;
        bool signFixed = false;
        real_t sign = -1;

        int_t minIndex = -1;
        int minZ = std::numeric_limits<int>::max(); // Initialize minZ to the maximum possible integer


        real_t kappa = 8.;
        real_t gamma = 5./3.;

        real_t m_i, Z;
        
        real_t PsiDiff(len_t);
        real_t InterpolatePsi(len_t);
        void FindRadiusOfLCFS();
        real_t StepFunction(len_t);
        void SetHeatLoss();
        void Deallocate();
        
	protected:
        /***************************** FUNCTIONS NEEDED FOR DIAGONALCOMPLEXTERM *****************************/
        virtual void SetWeights() override; 
        
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override; 
        /****************************************************************************************************/
        
        
    public:
    
        HaloRegionHeatLossTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*,
            real_t, bool userGivenPsiEdge_t0=0,
			real_t PsiEdge_t0=0
		);
        ~HaloRegionHeatLossTerm();
        // USES REBUILD, SETVECTORELEMENT, SETJACOBIANBLOCK FROM PARENT CLASSES
        
        virtual bool GridRebuilt() override;
        virtual const real_t *GetHaloRegionHeatLossWeights();
    
    };
}

#endif/*_DREAM_EQUATION_FLUID_PARALLEL_HEAT_LOSS_TERM_HPP*/
