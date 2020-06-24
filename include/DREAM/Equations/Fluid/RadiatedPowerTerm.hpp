#ifndef _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP
#define _DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP

#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "FVM/UnknownQuantityHandler.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/ADAS.hpp"

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
    class RadiatedPowerTerm : public FVM::DiagonalComplexTerm {
    private:
        ADAS *adas;
        IonHandler *ionHandler;

        len_t id_ncold;
        len_t id_Tcold;
        len_t id_ni;
    protected:
        virtual len_t GetNumberOfWeightsElements() override 
            {return ionHandler->GetNzs() * grid->GetNCells();}

        virtual void SetWeights() override;
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override;

    public:
        RadiatedPowerTerm(FVM::Grid*, FVM::UnknownQuantityHandler*, IonHandler*, ADAS*);
    };
}


#endif /*_DREAM_EQUATION_FLUID_RADIATED_POWER_TERM_HPP*/
