#ifndef _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP
#define _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"

namespace DREAM {
	class TrappingLimitedRRTransport :
	public RechesterRosenbluthTransport {
	protected:
		PitchScatterFrequency *nuD;
		len_t id_ncold, id_Tcold, id_ni;
		real_t *jacobianPrefactor = nullptr;
		real_t **drrDetrap = nullptr;
		real_t *kappaDetrap = nullptr;
		len_t np1_cache = 0;

		void AllocateCache();
		void ComputeKappa();

		virtual void SetPartialDiffusionTerm(len_t, len_t) override;

	public:
		TrappingLimitedRRTransport(
			FVM::Grid*, enum OptionConstants::momentumgrid_type,
			FVM::Interpolator1D*, CollisionQuantityHandler*,
				FVM::UnknownQuantityHandler*, bool withIonJacobian=true
		);
		virtual ~TrappingLimitedRRTransport();

		virtual bool GridRebuilt() override;
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		virtual bool SetJacobianBlock(const len_t, const len_t, FVM::Matrix*, const real_t*) override;
		const real_t *const* GetDetrappingLimit() const { return this->drrDetrap; }
		static real_t DetrappingEigenvalue(real_t xiT);
		const real_t *GetDetrappingEigenvalue() const { return this->kappaDetrap; }
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP*/
