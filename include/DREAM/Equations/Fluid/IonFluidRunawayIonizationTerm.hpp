#ifndef _DREAM_EQUATION_ION_FLUID_RUNAWAY_IONIZATION_TERM_HPP
#define _DREAM_EQUATION_ION_FLUID_RUNAWAY_IONIZATION_TERM_HPP


#include "DREAM/Constants.hpp"
#include "DREAM/ADAS.hpp"
#include "DREAM/Equations/Fluid/IonEquationTerm.hpp"
#include "DREAM/IonHandler.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/Equations/Fluid/IonKineticIonizationTerm.hpp"

namespace DREAM {
    class IonFluidRunawayIonizationTerm : public IonEquationTerm<FVM::EquationTerm> {
        private:

            len_t id_ions, id_nre;
            len_t tableIndexIon;

            FVM::UnknownQuantityHandler *u;
            real_t scaleFactor;

            real_t *weights;

            len_t nr;

            /**
             * The ionisation cross-section scales logarithmically with the runaway momentum p, 
             * and thus only weakly depends on p at relevant energies. This motivates the 
             * assumption of f_re ~ delta(p - p*) in approximating the runaway electron impact ionisation.
             * To this end, the characteristic momentum of 20 mc is used.
             */
            const real_t CHARACTERISTIC_RUNAWAY_MOMENTUM = 20.0;

	    real_t preFactor;




        public:
            IonFluidRunawayIonizationTerm(
                FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, const len_t, real_t
            );
            virtual ~IonFluidRunawayIonizationTerm();

            real_t *GetWeights() { return weights; }

            virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
            virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }


            virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;



            virtual bool SetCSJacobianBlock(
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

#endif/*_DREAM_EQUATION_ION_FLUID_RUNAWAY_IONIZATION_TERM_HPP*/
