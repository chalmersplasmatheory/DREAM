#ifndef _DREAM_EQUATIONS_FLUID_OHMIC_ELECTRIC_FIELD_TERM_HPP
#define _DREAM_EQUATIONS_FLUID_OHMIC_ELECTRIC_FIELD_TERM_HPP

#include "FVM/Equation/DiagonalLinearTerm.hpp"


namespace DREAM {
	class OhmicElectricFieldTerm : public FVM::DiagonalLinearTerm {
	protected:
		real_t scaleFactor;

		virtual void SetWeights() override;
	public:
		OhmicElectricFieldTerm(FVM::Grid*, const real_t);
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_OHMIC_ELECTRIC_FIELD_TERM_HPP*/
