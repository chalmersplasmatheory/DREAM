/**
 * Definition of equations relating to the poloidal flux.
 * In DREAM, the poloidal flux is normalised to the major
 * radius of the magnetic axis,
 * psi_p = poloidal flux / R0
 */

#include <iostream>
#include <string>
#include <gsl/gsl_sf_bessel.h>
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/Equation.hpp"
#include "FVM/Interpolator3D.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/AmperesLawDiffusionTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/HyperresistiveDiffusionTerm.hpp"

using namespace DREAM;
using namespace std;


/**
 * Implementation of a class which represents the j_||/(B/Bmin) term in Ampere's law.
 */
namespace DREAM {
    class AmperesLawJTotTerm : public FVM::DiagonalLinearTerm {
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
    public:
        AmperesLawJTotTerm(FVM::Grid* g) : FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = - Constants::mu0 * grid->GetRadialGrid()->GetFSA_1OverR2(ir) * grid->GetRadialGrid()->GetBTorG(ir) / grid->GetRadialGrid()->GetBmin(ir);
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}


#define MODULENAME "eqsys/psi_p"

/**
 * Construct the equation for the poloidal flux j_|| ~ Laplace(psi)
 * 
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquation_psi_p(
    EquationSystem *eqsys, Settings *s 
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Equation *eqn_j1 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn_j2 = new FVM::Equation(fluidGrid);


    // weightFunc3 represents -mu0*R0*<B*nabla phi>/Bmin (ie mu0*weightFunc1)
    std::function<real_t(len_t,len_t,len_t)> weightFunc3 = [fluidGrid](len_t ir,len_t, len_t)
        {return - Constants::mu0 * fluidGrid->GetRadialGrid()->GetFSA_1OverR2(ir) * fluidGrid->GetRadialGrid()->GetBTorG(ir) / fluidGrid->GetRadialGrid()->GetBmin(ir);};

    eqn_j1->AddTerm(new AmperesLawJTotTerm(fluidGrid));
    eqn_j2->AddTerm(new AmperesLawDiffusionTerm(fluidGrid));

    /**
     * TODO: Add additional boundary conditions.
     */

    eqsys->SetEquation(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_J_TOT, eqn_j1, "Poloidal flux Ampere's law");
    eqsys->SetEquation(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_POL_FLUX, eqn_j2);

    ConstructEquation_psi_p_initializeFromJ(eqsys,s);

}



/** 
 * Evaluates the initial psi_p(r) from an initial j_tot profile as the double integral
 *      psi(r) = psi(a) - 2*pi*mu0 * int( I_tot(r) / (VpVol*<|nabla r|^2/R^2>) , r, r_wall )
 * where
 *      I_tot = 1/(2*pi) * int( VpVol*j_tot(r) * G * <1/R^2>, 0, r )
 */

void SimulationGenerator::ConstructEquation_psi_p_initializeFromJ(EquationSystem *eqsys, Settings *s){
    FVM::RadialGrid *rGrid = eqsys->GetFluidGrid()->GetRadialGrid();
    len_t nr = rGrid->GetNr();
    real_t *psi_p_init = new real_t[nr];
    real_t *Itot = new real_t[nr];
    real_t *j_tot_init = LoadDataR("eqsys/j_tot", rGrid, s);

    const real_t *r = rGrid->GetR();
    const real_t *dr = rGrid->GetDr();
    #define integrand(I) 1/(2*M_PI) * rGrid->GetVpVol(I)*j_tot_init[I]*rGrid->GetBTorG(I)/rGrid->GetBmin(I) * rGrid->GetFSA_1OverR2(I)
    Itot[0] = r[0]*integrand(0);
    for(len_t ir=1; ir<nr; ir++){
        Itot[ir] = Itot[ir-1] + dr[ir-1]*integrand(ir);
    }
    #undef integrand

    const real_t rmax = rGrid->GetR_f(nr);
    #define integrand(I) 2*M_PI*Constants::mu0*Itot[I]/(rGrid->GetVpVol(I)*rGrid->GetFSA_NablaR2OverR2_f(I))
    psi_p_init[nr-1] = -(rmax-r[nr-1])*integrand(nr-1);
    for(len_t ir = nr-2; true; ir--){
        psi_p_init[ir] = psi_p_init[ir+1] - dr[ir]*integrand(ir);
        if(ir==0)
            break;
    }


}
