#ifndef _DREAM_EQUATIONS_KINETIC_TIME_VARYING_B_TERM_HPP
#define _DREAM_EQUATIONS_KINETIC_TIME_VARYING_B_TERM_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Interpolator1D.hpp"

namespace DREAM {
	class TimeVaryingBTerm : public FVM::AdvectionTerm {
	public:
		struct BAParams {
			real_t xi0;
		};
	private:
		real_t **BA_Fxi=nullptr;
		FVM::Interpolator1D *dlnB0dt;
		real_t current_dlnB0dt=0;

		static real_t BForce(
			real_t xiOverXi0, real_t BOverBmin, real_t ROverR0, 
			real_t NablaR2, void *par
		);

		void BounceAverageForce();

	public:
		TimeVaryingBTerm(
			FVM::Grid*, FVM::Interpolator1D*
		);
		~TimeVaryingBTerm();

		virtual bool GridRebuilt() override;
		virtual void Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler*) override;

		real_t **GetBounceAverage() { return this->BA_Fxi; }
		const real_t GetCurrentDbDt() const { return this->current_dlnB0dt; }
	};
}

#endif/*_DREAM_EQUATIONS_KINETIC_TIME_VARYING_B_TERM_HPP*/
