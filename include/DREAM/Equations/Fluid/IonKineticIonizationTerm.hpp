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

        OptionConstants::eqterm_ionization_mode ionization_mode;
        bool isPXiGrid;
        bool collfreqModeIsFull;
        len_t id_ions, id_nfast;
        real_t **IntegrandAllCS = nullptr;
        len_t tableIndexIon;
        real_t *tmpVec;
        
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
        void SetIntegrand(const len_t Z0, const len_t rOffset, real_t *diffIntegrand = nullptr);
        void RebuildIntegrand();
    public:
        IonKineticIonizationTerm(
            FVM::Grid*, FVM::Grid*, len_t momentId, len_t fId, FVM::UnknownQuantityHandler*, 
            IonHandler*, const len_t iIon, OptionConstants::eqterm_ionization_mode, 
            bool isPXiGrid, const len_t id_nf, 
            real_t pThreshold = 0, FVM::MomentQuantity::pThresholdMode pMode = FVM::MomentQuantity::P_THRESHOLD_MODE_MIN_MC
        );
        virtual ~IonKineticIonizationTerm();

        // nnz corresponding to one MomentQuantity per charge state (except the fully ionized one)
        virtual len_t GetNumberOfNonZerosPerRow() const override {
            return this->Zion * this->FVM::MomentQuantity::GetNumberOfNonZerosPerRow();
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

