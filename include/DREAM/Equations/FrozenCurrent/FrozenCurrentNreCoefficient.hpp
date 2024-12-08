#ifndef _DREAM_EQUATIONS_FROZEN_CURRENT_NRE_COEFFICIENT_HPP
#define _DREAM_EQUATIONS_FROZEN_CURRENT_NRE_COEFFICIENT_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTermHandler.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Interpolator1D.hpp"

namespace DREAM {
	class FrozenCurrentNreCoefficient : public FVM::EquationTerm {
	protected:
		FVM::Grid *fluidGrid = nullptr;
		FVM::Interpolator1D *I_p_presc = nullptr;
		RunawayFluid *REfluid = nullptr;

		real_t
			*gamma_gen = nullptr,
			*deriv = nullptr,
			*d2ndr2 = nullptr,
			**unit = nullptr;

		// Current value for D_I
		real_t D_I = 0;
		// Time for which the coefficient was rebuilt
		real_t coeffTime = -1;

		FVM::Operator
			*op_n_re,
			*op_n_tot,
			*op_n_i;

		FVM::DiffusionTerm *op_transport=nullptr;
		TransportDiffusiveBC *bc_transport=nullptr;

		len_t id_D_I, id_I_p, id_n_re, id_n_tot, id_n_i;
	
	public:
		FrozenCurrentNreCoefficient(
			FVM::Grid*, FVM::Grid*, FVM::Interpolator1D*,
			FVM::UnknownQuantityHandler*, RunawayFluid*,
			RunawaySourceTermHandler*
		);
		virtual ~FrozenCurrentNreCoefficient();

		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*);

		void SetTransportOperators(
			FVM::DiffusionTerm*, TransportDiffusiveBC*
		);

		virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
		virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
		virtual void SetVectorElements(real_t*, const real_t*) override;
	};
}

#endif/*_DREAM_EQUATIONS_FROZEN_CURRENT_NRE_COEFFICIENT_HPP*/
