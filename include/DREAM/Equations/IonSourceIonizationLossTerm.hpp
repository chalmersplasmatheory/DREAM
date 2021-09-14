#ifndef _DREAM_EQUATION_FLUID_ION_SOURCE_IONIZATION_LOSS_TERM_HPP
#define _DREAM_EQUATION_FLUID_ION_SOURCE_IONIZATION_LOSS_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NIST.hpp"

/**
 * Implementation of an equation term describing the energy loss required to ionize 
 * an ion source to the charge state distribution it has when added to the plasm
 */
namespace DREAM {
    class IonSourceIonizationLossTerm : public FVM::EquationTerm {
    private:
        EquationTerm *ionSource;
        IonHandler *ionHandler;
        NIST *nist;
        static const len_t nZ;
        static const len_t nData;
        static const real_t dataZ[];
        static const real_t dataIbind[];

        real_t *ETotIoniz;

        void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler);

    public:
        IonSourceIonizationLossTerm(FVM::Grid* g, IonHandler *ionHandler, NIST *nist, EquationTerm *ionSource);

        virtual void SetMatrixElements(FVM::Matrix* mat, real_t*) override; 
        virtual void SetJacobianBlock(FVM::Matrix* jac, real_t*) override;
        virtual void SetVectorElements(real_t* vec, const real_t*) override;
    };
}


#endif /*_DREAM_EQUATION_FLUID_ION_SOURCE_IONIZATION_LOSS_TERM_HPP*/
