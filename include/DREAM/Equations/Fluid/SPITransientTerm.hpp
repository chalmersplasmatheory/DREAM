#ifndef _DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP
#define _DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP

#include "FVM/config.h"
#include "FVM/Equation/DiagonalTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Matrix.hpp"


namespace DREAM {
    class SPITransientTerm : public FVM::DiagonalTerm {
    private:
        real_t dt;

        // ID of differentiated quantity
        len_t unknownId;
        len_t nShard;
        real_t scaleFactor;

        // Differentiated quantity at the previous time step
        real_t *xn;
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
        virtual void AddWeightsJacobian(const len_t, const len_t, FVM::Matrix*, const real_t*) override {}

        virtual len_t GetNumberOfWeightsElements() override
            {return this->nShard;}
    public:
        SPITransientTerm(FVM::Grid*, const len_t unknownId, len_t nShard, real_t scaleFactor=1.0);
        virtual void SetWeights() override;
        virtual void Rebuild(const real_t, const real_t dt, FVM::UnknownQuantityHandler *uqty) override;
        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
    };
}

#endif/*_DREAM_EQUATION_FLUID_SPI_TRANSIENT_TERM_HPP*/
