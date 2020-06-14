#ifndef _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP
#define _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"

/**
 * Implementation of a class which represents the 
 * radiated power as calculated with rate coefficients
 * from the ADAS database (PLT corresponds to line
 * and PRB to brems and recombination radiated power).
 * The term is of the form n_e * sum_i n_i L_i, summed over all
 * ion species i. In the semi-implicit solver, n_e is the "unknown"
 * evaluated at the next time step and n_i L_i coefficients.
 * We ignore the Jacobian with respect to L_i(n,T) and capture only the
 * n_e and n_i contributions.
 */
namespace DREAM {
    class RadiatedPowerTerm : public FVM::DiagonalQuadraticTerm {
    private:
        ADAS *adas;
        IonHandler *ionHandler;
    protected:
        // radiated power coefficients depends on T (and weakly on n)
        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void SetWeights() override
            {
                len_t NCells = grid->GetNCells();
                len_t nZ = ionHandler->GetNZ();
                const len_t *Zs = ionHandler->GetZs();

                len_t id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
                len_t id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
                real_t *n_cold = unknowns->GetUnknownData(id_ncold);
                real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
                for(len_t iz = 0; iz<nZ; iz++){
                    ADASRateInterpolator *PLT_interper = adas->GetPLT(Zs[iz]);
                    ADASRateInterpolator *PRB_interper = adas->GetPRB(Zs[iz]);
                    for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){

                        len_t offset = 0;
                        len_t nMultiple = ionHandler->GetIndex(iz,Z0);
                        for (len_t ir = 0; ir < nr; ir++){
                            real_t w =  PLT_interper->Eval(Z0, n_cold[ir], T_cold[ir])
                                      + PRB_interper->Eval(Z0, n_cold[ir], T_cold[ir]);
                            for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                                    weights[NCells*nMultiple + offset + i] = w;
                            offset += n1[ir]*n2[ir];
                        }
                    }
                }
            }
    public:
        RadiatedPowerTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, IonHandler *ionHandler, ADAS *adas) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES),u) 
        {
            this->adas = adas;
            this->ionHandler = ionHandler;
        }
    };
}


#endif /*_DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP*/
