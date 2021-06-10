#ifndef _DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP
#define _DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM {
    class SPITransientTerm : public FVM::EquationTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        len_t nShard;
        real_t scaleFactor;

        // Differentiated quantity at the previous time step
        real_t *xn;
    public:
        SPITransientTerm(FVM::Grid*, const len_t unknownId, len_t nShard, real_t scaleFactor=1.0);
        ~SPITransientTerm(){}

        virtual len_t GetNumberOfNonZerosPerRow() const { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const { return 1; }

        virtual void Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler *uqty) override;
        virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP*/
