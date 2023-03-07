#ifndef _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_HPP
#define _DREAM_EQUATIONS_FROZEN_CURRENT_COEFFICIENT_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Interpolator1D.hpp"

namespace DREAM {
	class FrozenCurrentCoefficient : public FVM::EquationTerm {
	protected:
		FVM::Grid *fluidGrid = nullptr;
		FVM::Interpolator1D *I_p_presc=nullptr;
		len_t id_D_I, id_I_p, id_j_tot;
		real_t prevTime=-1;

		const int nMaxPoints = 4;
		int nPoints = 0;
		// Estimates of diffusion coefficient
		real_t *D_I, DI_sol;
		//real_t D_I = 0, D_I_prev = 0;
		real_t D_I_max = 1e3;

		// Start by not prescribing plasma current
		real_t Ip_presc;
		real_t *Ip;

		real_t EstimateSolutionZeroth(FVM::UnknownQuantityHandler*, const real_t, const real_t);
		real_t EstimateSolutionLinear(const real_t);
		real_t EstimateSolutionQuadratic(const real_t);

		void ClearSolutions();
		void PushSolution(const real_t, const real_t);

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
