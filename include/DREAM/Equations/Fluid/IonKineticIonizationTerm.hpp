#ifndef _DREAM_EQUATION_ION_KINETIC_IONIZATION_TERM_HPP
#define _DREAM_EQUATION_ION_KINETIC_IONIZATION_TERM_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/MomentQuantity.hpp"

namespace DREAM {
    class IonKineticIonizationTerm : public IonEquationTerm<FVM::MomentQuantity> {
    private:
        struct kinetic_ionization_rate {
            const char *name;
            len_t Z;
            const real_t *params; // parameters in the model (size nParams x Z)
        };

        len_t id_ions;
        bool isPXiGrid;
        real_t **IntegrandAllCS = nullptr;
        len_t tableIndexIon;
        
        // number of parameters used in the fit
        static const len_t nParamsForFit;
        // ionization rate atomic data
        static const len_t kinetic_rate_n;
        static struct kinetic_ionization_rate kinetic_rate_table[];


        static real_t EvaluateBCGSingleSubshell(real_t p, real_t C, real_t I_pot_eV, real_t betaStar);
        static real_t EvaluateIonizationCrossSection(real_t p, const real_t *params);
        static int_t GetTableIndex(len_t Z);

        void Allocate();
        void Deallocate();
        void SetIntegrand(const len_t Z0, const len_t rOffset, real_t *diffIntegrand);
        void RebuildIntegrand();
    public:
        IonKineticIonizationTerm(
            FVM::Grid*, FVM::Grid*, len_t momentId, len_t fId, FVM::UnknownQuantityHandler*, 
            IonHandler*, const len_t iIon, bool isPXiGrid
        );
        virtual ~IonKineticIonizationTerm();

        //static real_t EvaluateIonizationCrossSection(real_t p, const len_t Z, const len_t Z0);


        //virtual len_t GetNumberOfNonZerosPerRow() const override 
        //    { return this->FVM::MomentQuantity::GetNumberOfNonZerosPerRow(); }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            {
                len_t nnz = this->GetNumberOfNonZerosPerRow();
                nnz += 2; // 1 for ncold partial derivative and 1 for Tcold 
                return nnz; 
            }

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override{}

        virtual void SetCSJacobianBlock(
            const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t *nions,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSMatrixElements(
            FVM::Matrix *mat, real_t *rhs, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
        ) override;
        virtual void SetCSVectorElements(
            real_t *vec, const real_t *nions, const len_t /*iIon*/, const len_t Z0, const len_t rOffset
        ) override;
    };
}

#endif/*_DREAM_EQUATION_ION_KINETIC_IONIZATION_TERM_HPP*/

