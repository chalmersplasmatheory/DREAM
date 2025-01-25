#ifndef _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP
#define _DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class HeatTransportRechesterRosenbluth : public FVM::DiffusionTerm {
    protected:
        FVM::Interpolator1D *deltaBOverB=nullptr;
        FVM::UnknownQuantityHandler *unknowns;
        real_t *dD=nullptr;

        // IDs of unknown quantities used by the operator...
        len_t id_n_cold, id_T_cold;

        void AllocateDiffCoeff();
        virtual void SetPartialDiffusionTerm(len_t, len_t) override;

		virtual const real_t *EvaluateDeltaBOverB(const real_t t) { return this->deltaBOverB->Eval(t); }

    public:
        HeatTransportRechesterRosenbluth(FVM::Grid*, FVM::Interpolator1D*, FVM::UnknownQuantityHandler*);
        virtual ~HeatTransportRechesterRosenbluth();

        virtual bool GridRebuilt() override;
        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_FLUID_HEAT_TRANSPORT_RECHESTER_ROSENBLUTH_HPP*/
