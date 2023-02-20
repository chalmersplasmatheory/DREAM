#ifndef _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_I_HPP
#define _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_I_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Interpolator1D.hpp"

namespace DREAM {
	class FrozenCurrentCoefficient_I : public FVM::EquationTerm {
	protected:
		FVM::Grid *fluidGrid = nullptr;
		FVM::Interpolator1D *I_p_presc=nullptr;
		len_t id_I_p;
		real_t prevTime=-1;

		real_t Ip, Ipresc;
		bool prescribeCurrent = false;

	public:
		FrozenCurrentCoefficient_I(
			FVM::Grid*, FVM::Interpolator1D*,
			FVM::UnknownQuantityHandler*
		);
		virtual ~FrozenCurrentCoefficient_I();

		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

		virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;
	};
}

#endif/*_DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_I_HPP*/
