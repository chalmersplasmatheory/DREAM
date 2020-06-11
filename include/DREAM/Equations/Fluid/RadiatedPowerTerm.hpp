#ifndef _DREAM_EQUATION_RADIATED_POWER_TERM_HPP
#define _DREAM_EQUATION_RADIATED_POWER_TERM_HPP

#include "FVM/Equation/DiagonalQuadraticTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"

/**
 * Implementation of a class which represents the 
 * radiated power as calculated with rate coefficients
 * from the ADAS database (retrieved by GetRadiatedPowerCoefficient).
 * Term is of the form n_e * sum_i n_i L_i, summed over all
 * ion species i. In the semi-implicit solver, n_e is the "unknown"
 * evaluated at the next time step and n_i L_i coefficients.
 * We ignore the Jacobian with respect to L_i(n,T) and capture only the
 * n_e and n_i contributions.
 */
namespace DREAM {
    class RadiatedPowerTerm : public FVM::DiagonalQuadraticTerm {
    private:
        IonHandler *ionHandler;
        real_t GetRadiatedPowerCoefficient(len_t, len_t, real_t, real_t)
            {return 0;}
    protected:
        // radiated power coefficients depend on T (and weakly on n)
        virtual bool TermDependsOnUnknowns() override {return true;}
    public:
        RadiatedPowerTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, IonHandler *ionHandler) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES),u), ionHandler(ionHandler){}

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
                for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){

                    len_t offset = 0;
                    len_t nMultiple = ionHandler->GetIndex(iz,Z0);
                    for (len_t ir = 0; ir < nr; ir++){
                        real_t w = GetRadiatedPowerCoefficient(iz,Z0,n_cold[ir],T_cold[ir]);
                        for(len_t i = 0; i < n1[ir]; i++)
                            for(len_t j = 0; j < n2[ir]; j++)
                                weights[NCells*nMultiple + offset + n1[ir]*j + i] = w;
                        offset += n1[ir]*n2[ir];
                    }
                }
            }
        }
    };
}


#endif /*_DREAM_EQUATION_RADIATED_POWER_TERM_HPP*/
