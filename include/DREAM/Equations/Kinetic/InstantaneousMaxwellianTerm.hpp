#ifndef _DREAM_EQUATIONS_KINETIC_INSTANTANEOUS_MAXWELLIAN_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_INSTANTANEOUS_MAXWELLIAN_TERM_HPP

#include "FVM/Equation/EquationTerm.hpp"


namespace DREAM {
	class InstantaneousMaxwellianTerm : public FVM::EquationTerm {
	protected:
		len_t id_T, id_n;
		real_t *F, *dFdn, *dFdT;

	public:
		enum MaxwellianPopulation {
			MAXWELLIAN_POPULATION_COLD,
			MAXWELLIAN_POPULATION_HOT
		};

		InstantaneousMaxwellianTerm(
			FVM::Grid*, FVM::UnknownQuantityHandler*,
			enum MaxwellianPopulation
		);
		virtual ~InstantaneousMaxwellianTerm();

        virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
        virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return GetNumberOfNonZerosPerRow(); }

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

		virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_INSTANTANEOUS_MAXWELLIAN_TERM_HPP*/
