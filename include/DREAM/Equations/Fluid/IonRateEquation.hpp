#ifndef _DREAM_EQUATION_ION_RATE_EQUATION_HPP
#define _DREAM_EQUATION_ION_RATE_EQUATION_HPP

#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class IonRateEquation : public IonEquationTerm<FVM::EquationTerm> {
    private:
        enum SetMode {MATRIX, JACOBIAN};
        ADAS *adas;
        FVM::UnknownQuantityHandler *unknowns;
        len_t id_ions, id_n_cold, id_n_tot, id_T_cold;
        bool addFluidIonization; // the full ADAS ionization rate is added in this equation term
        bool addFluidJacobian;   // only the jacobian of the ionization is set with this term
        real_t
            **Rec,         // Radiative recombination rates (nZs x nr)
            **PartialNRec, // d/dn_cold of radiative recombination rates  (nZs x nr)
            **PartialTRec, // d/dT_cold of radiative recombination rates  (nZs x nr)
            **Ion,         // Ionization rate coefficients (nZs x nr)
            **PartialNIon, // d/dn_cold of ionization rate coefficients (nZs x nr)
            **PartialTIon; // d/dT_cold of ionization rate coefficients (nZs x nr)
    public:
        IonRateEquation(
            FVM::Grid*, IonHandler*, const len_t, ADAS*, 
            FVM::UnknownQuantityHandler*,bool,bool
        );
        virtual ~IonRateEquation();

        void AllocateRateCoefficients();
        void DeallocateRateCoefficients();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 3; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override 
            {
                len_t nnz = this->GetNumberOfNonZerosPerRow();
                nnz += 2; // 1 for ncold partial derivative and 1 for Tcold 
                return nnz; 
            }

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

        virtual bool SetCSJacobianBlock(
            const len_t, const len_t, FVM::Matrix*, const real_t*,
            const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;

        virtual void SetCSMatrixElements(
            FVM::Matrix *mat, real_t *rhs, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override {
            SetCSMatrixElements(mat, rhs, iIon, Z0, rOffset, MATRIX);
        }

        virtual void SetCSMatrixElements(
            FVM::Matrix*, real_t*, const len_t iIon, const len_t Z0, const len_t rOffset, SetMode
        );

        virtual void SetCSVectorElements(
            real_t*, const real_t*, const len_t iIon, const len_t Z0, const len_t rOffset
        ) override;
    };
}

#endif/*_DREAM_EQUATION_ION_RATE_EQUATION_HPP*/

