#ifndef _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP
#define _DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP

#include "DREAM/Equations/Kinetic/RechesterRosenbluthTransport.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

namespace DREAM {
	class TrappingLimitedRRTransport :
		public RechesterRosenbluthTransport {
	protected:
		RunawayFluid *REFluid;

	public:
		TrappingLimitedRRTransport(
			FVM::Grid*, enum OptionConstants::momentumgrid_type,
			FVM::Interpolator1D*, RunawayFluid*
		);
		virtual ~TrappingLimitedRRTransport();

		virtual void Rebuild() override;
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_TRAPPING_LIMITED_RR_TRANSPORT_HPP*/
