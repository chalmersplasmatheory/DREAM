#ifndef _DREAM_EQUATION_FLUID_NET_ION_DENSITY_FROM_ION_CHARGE_STATES_TERM_HPP
#define _DREAM_EQUATION_FLUID_NET_ION_DENSITY_FROM_ION_CHARGE_STATES_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"

namespace DREAM {
    class NetIonDensityFromIonChargeStatesTerm : public FVM::EquationTerm {
    private:
        const len_t Z;
        const len_t iz;
        IonHandler *ionHandler;
        real_t scaleFactor;
    public:
        NetIonDensityFromIonChargeStatesTerm(FVM::Grid *g, const len_t Z, const len_t iz, IonHandler *ionHandler, real_t scaleFactor=1.0) 
            : FVM::EquationTerm(g), Z(Z), iz(iz), ionHandler(ionHandler), scaleFactor(scaleFactor){}
        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override {}

        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t*) override {
            if(uqtyId==derivId)
                SetMatrixElements(jac, nullptr);
        }
        virtual void SetMatrixElements(FVM::Matrix *mat, real_t *) override {
            for(len_t Z0=0; Z0<=Z; Z0++){
                len_t indZ = ionHandler->GetIndex(iz,Z0);
                for(len_t ir=0; ir<nr; ir++)
                    mat->SetElement(iz*nr+ir, indZ*nr+ir, scaleFactor);
            }
        }
        virtual void SetVectorElements(real_t *vec, const real_t *ni) override {
            for(len_t Z0=0; Z0<=Z; Z0++){
                len_t indZ = ionHandler->GetIndex(iz,Z0);
                for(len_t ir=0; ir<nr; ir++)
                    vec[iz*nr+ir] += scaleFactor*ni[indZ*nr+ir];
            }
        }
    };
}

#endif/*_DREAM_EQUATION_FLUID_NET_ION_DENSITY_FROM_ION_CHARGE_STATES_TERM_HPP*/
