#ifndef _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP
#define _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP

#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"

namespace DREAM {
	class TrappingLimitedRRTransport :
	public RechesterRosenbluthTransport {
	protected:
		PitchScatterFrequency *nuD;
		enum OptionConstants::eqterm_transport_rechester_rosenbluth_detrapping_mode detrappingMode;
		real_t detrappingScale;

		len_t id_ncold, id_Tcold, id_ni;
		real_t *jacobianPrefactor = nullptr;
		real_t **drrDetrap = nullptr;
		len_t np1_cache = 0;

		void AllocateCache();

		virtual void SetPartialDiffusionTerm(len_t, len_t) override;

	public:
		TrappingLimitedRRTransport(
			FVM::Grid*, enum OptionConstants::momentumgrid_type,
			FVM::Interpolator1D*, CollisionQuantityHandler*,
			FVM::UnknownQuantityHandler*, bool withIonJacobian=true,
			enum OptionConstants::eqterm_transport_rechester_rosenbluth_detrapping_mode
				detrappingMode=OptionConstants::EQTERM_TRANSPORT_RECHESTER_ROSENBLUTH_DETRAPPING_MODE_LOCAL,
			real_t detrappingScale=1.0
		);
		virtual ~TrappingLimitedRRTransport();

		virtual bool GridRebuilt() override;
		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
		const real_t *const* GetDetrappingLimit() const { return this->drrDetrap; }
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP*/
