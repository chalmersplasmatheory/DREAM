#ifndef _DREAM_EQUATIONS_FROZEN_CURRENT_NRE_COEFFICIENT_HPP
#define _DREAM_EQUATIONS_FROZEN_CURRENT_NRE_COEFFICIENT_HPP

#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/Equations/RunawaySourceTermHandler.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/Matrix.hpp"

namespace DREAM {
	class FrozenCurrentNreCoefficient : public FVM::EquationTerm {
	protected:
		FVM::Grid *fluidGrid = nullptr;
		FVM::Interpolator1D *I_p_presc = nullptr;
		RunawayFluid *REfluid = nullptr;
		FVM::UnknownQuantityHandler *unknowns = nullptr;

		real_t t_adjust = 0.01;

		real_t
			S_gen, S_loss,
			*gamma_gen = nullptr,
			*d2ndr2 = nullptr,
			*dx_vec = nullptr,
			*dSgen_dx = nullptr,
			*dSloss_dx = nullptr,
			**unit = nullptr;

		FVM::Matrix *op_mat;
		len_t nMatrixElements;

		// Current value for D_I
		real_t D0, dD;
		real_t D_I = 0;
		real_t lastDt, dIp, dIohm_dt;

		FVM::Operator
			*op_n_re,
			*op_n_tot,
			*op_n_i;

		FVM::DiffusionTerm *op_transport=nullptr;
		TransportDiffusiveBC *bc_transport=nullptr;

		len_t id_D_I, id_I_p, id_n_re, id_n_tot, id_n_i, id_j_ohm;

		real_t EvaluateCurrentIntegral(const real_t*, bool densityToCurrent=true);
	
	public:
		FrozenCurrentNreCoefficient(
			FVM::Grid*, FVM::Grid*, FVM::Interpolator1D*,
			FVM::UnknownQuantityHandler*, RunawayFluid*,
			RunawaySourceTermHandler*, const real_t
		);
		virtual ~FrozenCurrentNreCoefficient();

		virtual len_t GetNumberOfNonZerosPerRow() const override { return 1; }
		virtual len_t GetNumberOfNonZerosPerRow_jac() const override { return 1; }

		real_t GetD0() { return this->D0; }
		real_t GetDD() { return this->dD; }
		real_t GetDIp() { return this->dIp; }
		real_t GetSgen() { return this->S_gen; }
		real_t GetSloss() { return this->S_loss; }

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
