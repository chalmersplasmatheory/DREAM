#ifndef _DREAM_EQUATIONS_FLUID_RUNAWAY_TRANSPORT_RECHESTER_ROSENBLUTH_HPP
#define _DREAM_EQUATIONS_FLUID_RUNAWAY_TRANSPORT_RECHESTER_ROSENBLUTH_HPP

#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
	class RunawayTransportRechesterRosenbluth : public FVM::DiffusionTerm {
	protected:
		FVM::Interpolator1D *dBB;

		virtual const real_t *EvaluateDeltaBOverB(const real_t t) { return dBB->Eval(t); }

	public:
		RunawayTransportRechesterRosenbluth(
			FVM::Grid*, FVM::Interpolator1D*
		);
		virtual ~RunawayTransportRechesterRosenbluth();

		virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*);
	};
}

#endif/*_DREAM_EQUATIONS_FLUID_RUNAWAY_TRANSPORT_RECHESTER_ROSENBLUTH_HPP*/
