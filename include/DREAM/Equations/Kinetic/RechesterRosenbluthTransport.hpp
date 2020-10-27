#ifndef _DREAM_EQUATIONS_RECHESTER_ROSENBLUTH_HPP
#define _DREAM_EQUATIONS_RECHESTER_ROSENBLUTH_HPP

#include "DREAM/Settings/OptionConstants.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "FVM/Interpolator1D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    class RechesterRosenbluthTransport : public FVM::DiffusionTerm {
    private:
        enum OptionConstants::momentumgrid_type mgtype;
        FVM::Interpolator1D *deltaBOverB;

        FVM::UnknownQuantityHandler *unknowns;

    public:
        RechesterRosenbluthTransport(FVM::Grid*, enum OptionConstants::momentumgrid_type, FVM::Interpolator1D*);
        ~RechesterRosenbluthTransport();

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
    };
}

#endif/*_DREAM_EQUATIONS_RECHESTER_ROSENBLUTH_HPP*/
