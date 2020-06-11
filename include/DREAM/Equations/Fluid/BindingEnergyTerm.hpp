#ifndef _DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP
#define _DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "DREAM/IonHandler.hpp"

/**
 * Implementation of a class which represents total potential (binding) 
 * energy of the ions in the plasma.
 */
namespace DREAM {
    class BindingEnergyTerm : public FVM::DiagonalLinearTerm {
    private:
        IonHandler *ionHandler;

        /**
         * TODO: Implement this method which returns the binding energy
         * of ion with atomic number Z and charge Z0.
         */
        real_t GetBindingEnergy(len_t /*Z*/, len_t /*Z0*/){
            return 0;
        }
    protected:
        virtual len_t GetNumberOfWeightsElements() override
            {return ionHandler->GetNzs() * grid->GetNCells();}


        virtual void SetWeights() override {
            len_t N = grid->GetNCells();
            len_t nZ = ionHandler->GetNZ();
            const len_t *Zs = ionHandler->GetZs();
            for(len_t iz = 0; iz<nZ; iz++){
                for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                    len_t n = ionHandler->GetIndex(iz,Z0);
                    real_t w = GetBindingEnergy(Zs[iz],Z0);
                    for(len_t i = 0; i<N; i++)
                        weights[n*N+i] = w;
                }
            }
        }
    public:
        BindingEnergyTerm(FVM::Grid* g, IonHandler *ionHandler) 
            : FVM::DiagonalLinearTerm(g), ionHandler(ionHandler){}

        virtual void SetMatrixElements(FVM::Matrix* mat, real_t*) override 
            { 
                len_t N = grid->GetNCells();
                len_t nnz = ionHandler->GetNzs();
                for(len_t i; i<N*nnz; i++)
                    mat->SetElement(i,i,weights[i]);

            }
        virtual void SetVectorElements(real_t* vec, const real_t* x) override 
            {
                len_t N = grid->GetNCells();
                len_t nnz = ionHandler->GetNzs();
                for(len_t i; i<N*nnz; i++)
                    vec[i] += weights[i] * x[i];
            }

    };
}


#endif /*_DREAM_EQUATION_FLUID_BINDING_ENERGY_TERM_HPP*/
