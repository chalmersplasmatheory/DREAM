#ifndef _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_HPP
#define _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Interpolator1D.hpp"

namespace DREAM {
	class FrozenCurrentCoefficient : public FVM::EquationTerm {
	protected:
		FVM::Grid *fluidGrid = nullptr;
		FVM::Interpolator1D *I_p_presc=nullptr;
		len_t id_I_p, id_j_tot;
		real_t prevTime=-1;

		// Current estimate of diffusion coefficient
		real_t D_I = 0, D_I_prev = 0;
		real_t D_I_max = 1e3;

		// Start by not prescribing plasma current
		real_t Ip, Ip_presc, Ip_prev;

	public:
		FrozenCurrentCoefficient(
			FVM::Grid*, FVM::Grid*, FVM::Interpolator1D*,
			FVM::UnknownQuantityHandler*, const real_t D_I_max=1e3
		);
		virtual ~FrozenCurrentCoefficient();

		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;

		virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;
	};
}

#endif/*_DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_HPP*/
