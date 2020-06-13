#ifndef _DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP
#define _DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NIST.hpp"

/**
 * Implementation of a class which represents total potential (binding) 
 * energy of the ions in the plasma.
 */
namespace DREAM {
    class BindingEnergyTerm : public FVM::DiagonalLinearTerm {
    private:
        IonHandler *ionHandler;
        NIST *nist;
        static const len_t nZ;
        static const len_t nData;
        static const real_t dataZ[];
        static const real_t dataIbind[];

        /*real_t GetBindingEnergy(len_t Zin, len_t Z0in);
        bool DataSetIsValid();*/
    protected:
        virtual len_t GetNumberOfWeightsElements() override
            {return ionHandler->GetNzs() * grid->GetNCells();}
        virtual void SetWeights() override;

    public:
        BindingEnergyTerm(FVM::Grid* g, IonHandler *ionHandler, NIST *nist);

        virtual void SetMatrixElements(FVM::Matrix* mat, real_t*) override; 
        virtual void SetVectorElements(real_t* vec, const real_t* x) override;
    };
}


#endif /*_DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP*/
