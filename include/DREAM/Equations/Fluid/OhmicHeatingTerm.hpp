#ifndef _DREAM_EQUATION_FLUID_OHMIC_HEATING_TERM_HPP
#define _DREAM_EQUATION_FLUID_OHMIC_HEATING_TERM_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"

/**
 * Implementation of a class which represents the 
 * ohmic heating term <j_ohm*E_||>.
 * In the semi-implicit solver we view j_ohm as the fixed 
 * (coefficient) quantity, and Eterm as the "unknown".
 */
namespace DREAM {
    class OhmicHeatingTerm : public FVM::DiagonalQuadraticTerm {
    protected:
        virtual void SetWeights() override
        {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                        weights[offset + i] = -w;
                offset += n1[ir]*n2[ir];
            }
        }

    public:
        OhmicHeatingTerm(
			FVM::Grid* g, FVM::UnknownQuantityHandler *u,
			const char *uqty_j=OptionConstants::UQTY_J_OHM
		) : FVM::DiagonalQuadraticTerm(g, u->GetUnknownID(uqty_j), u){}

		OhmicHeatingTerm(
			FVM::Grid* g, FVM::UnknownQuantityHandler *u,
			const len_t id_j
		) : FVM::DiagonalQuadraticTerm(g, id_j, u) {}
    };
}


#endif /*_DREAM_EQUATION_FLUID_OHMIC_HEATING_TERM_HPP*/
