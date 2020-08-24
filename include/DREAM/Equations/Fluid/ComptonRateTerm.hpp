#ifndef _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP
#define _DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "FVM/Equation/DiagonalComplexTerm.hpp"
//#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
/**
 * Implementation of a class which represents the Gamma_tritium contribution to the n_re equation.
 * Employs the analytical growth rate calculated by RunawayFluid.
 */
namespace DREAM {
    class ComptonRateTerm : public FVM::EquationTerm {
    private:
        RunawayFluid *REFluid;
        len_t id_n_tot;
        real_t scaleFactor;
        real_t *gamma;
        const real_t *data_n_tot;
    protected: 
    public:
        ComptonRateTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
                RunawayFluid *ref, real_t scaleFactor = 1.0);
        ~ComptonRateTerm();

        void AllocateGamma();
        void DeallocateGamma();
        
        virtual bool GridRebuilt() override;
        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 1; }   /* XXX TODO XXX */
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual void SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;

    };
}

#endif /*_DREAM_EQUATION_FLUID_COMPTON_RATE_TERM_HPP*/
