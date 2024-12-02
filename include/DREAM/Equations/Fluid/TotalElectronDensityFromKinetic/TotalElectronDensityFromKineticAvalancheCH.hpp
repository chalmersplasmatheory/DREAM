#ifndef _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_CH_HPP
#define _DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_CH_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"

/**
 * Implementation of an equation term which represents the total
 * number of electrons created by the kinetic Chiu-Harvey source
 */ 
namespace DREAM {
    class TotalElectronDensityFromKineticAvalancheCH 
        : public FVM::DiagonalComplexTerm {
    public:
        real_t pCutoff;
        AvalancheSourceCH *avaCH;
        real_t scaleFactor;
        len_t id_fhot;
        TotalElectronDensityFromKineticAvalancheCH(FVM::Grid* grid, FVM::Grid* og, real_t pCutoff, FVM::UnknownQuantityHandler *u, AvalancheSourceCH *avaCH, real_t scaleFactor = 1.0);

        virtual void SetWeights() override;
        
        virtual void SetDiffWeights(len_t, len_t) override {}

        bool AddWeightsJacobian(
            const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t* x
        ) override;

        /**
        * Set the linear operator matrix elements corresponding to this term.
        */
        void SetMatrixElements(FVM::Matrix *mat, real_t*) override;

        /**
        * Set function vector for this term.
        */
        void SetVectorElements(real_t *vec, const real_t *x) override;


        bool SetJacobianBlock(
            const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x
        )  override;
    };
}

#endif/*_DREAM_EQUATIONS_TOTAL_ELECTRON_DENSITY_FROM_KINETIC_AVALANCHE_SOURCE_CH_HPP*/