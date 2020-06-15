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
#include "FVM/Equation/Operator.hpp"
#include "FVM/Interpolator3D.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/AmperesLawDiffusionTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/HyperresistiveDiffusionTerm.hpp"
#include "DREAM/Equations/Scalar/WallCurrentTerms.hpp"

using namespace DREAM;
using namespace std;


/**
 * Implementation of a class which represents the j_||/(B/Bmin) term in Ampere's law.
 */
namespace DREAM {
    class AmperesLawJTotTerm : public FVM::DiagonalLinearTerm {
    public:
        AmperesLawJTotTerm(FVM::Grid* g) : FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = - Constants::mu0 * grid->GetRadialGrid()->GetFSA_1OverR2(ir) * grid->GetRadialGrid()->GetBTorG(ir) / grid->GetRadialGrid()->GetBmin(ir);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
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
    eqsys->SetUnknown(OptionConstants::UQTY_POL_FLUX, fluidGrid);

    FVM::Operator *eqn_j1 = new FVM::Operator(fluidGrid);
    FVM::Operator *eqn_j2 = new FVM::Operator(fluidGrid);

    eqn_j1->AddTerm(new AmperesLawJTotTerm(fluidGrid));
    eqn_j2->AddTerm(new AmperesLawDiffusionTerm(fluidGrid));

    /**
     * TODO: Add additional boundary conditions.
     */

    eqsys->SetOperator(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_J_TOT, eqn_j1, "Poloidal flux Ampere's law");
    eqsys->SetOperator(OptionConstants::UQTY_POL_FLUX, OptionConstants::UQTY_POL_FLUX, eqn_j2);

    /**
     * Initialization: define the function which integrates j_tot.
     * In principle, this expression is a direct inversion of the 
     * equation given previously with boundary condition psi(r_max) = 0 
     * It would be nicer to just solve the equation. (may be tricky since 
     * the boundary condittion may differ from the one we wish to use later)
     */
    const len_t id_psi = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_POL_FLUX);
    const len_t id_jtot = eqsys->GetUnknownHandler()->GetUnknownID(OptionConstants::UQTY_J_TOT);
    FVM::RadialGrid *rGrid = fluidGrid->GetRadialGrid();
    std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_PsiPFromJtot 
        = [rGrid,id_jtot](FVM::UnknownQuantityHandler*u, real_t *psi_p_init)
        {
            len_t nr = rGrid->GetNr();
            real_t *Itot = new real_t[nr];

            real_t *j_tot_init = u->GetUnknownData(id_jtot);
            
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
            #undef integrand
            delete [] Itot;
        };

    eqsys->initializer->AddRule(
        id_psi,
        EqsysInitializer::INITRULE_EVAL_FUNCTION,
        initfunc_PsiPFromJtot,
        // Dependencies
        id_jtot
    );

    ConstructEquations_I_wall(eqsys,s);


}

/**
 * Construct the equations for total toroidal plasma current, 
 * wall current and wall poloidal flux. Since poloidal flux is
 * arbitrary up to a constant additative factor, we choose the 
 * convention that the initial wall flu psi_w is always = 0.
 * The three coupled scalar equations to solve are:
 *   eq1: 0 = dpsi_w/dt + R_w I_w
 *   eq2: 0 = psi_w - L_w*I_w - psi(r_max) - I_p * integral_term
 *   eq3: 0 = I_p - integral(j_tot)
 *
 * 
 * eqsys: Equation system to put the equation in.
 * s:     Settings object describing how to construct the equation.
 */
void SimulationGenerator::ConstructEquations_I_wall(
    EquationSystem *eqsys, Settings */*s*/ 
) {
    // From settings, get wall resistivity R_W and wall inductance L_W.
    // Typical case is R_W = 0, in which case it is perfectly conducting 
    // and dpsi/dt = Vloop = 0 on the wall, which corresponds to the GO case.

    // Perhaps fix L_W to a characteristic inductance and use tau_wall as the only input?
    real_t L_W = 1e-6; // ~mu0 (inductance normalized to R0)
    real_t R_W = L_W * 100; // L_W * 1/tau_wall, tau_wall ~ 10ms in DIII-D.

    real_t a = 1; //Placeholder: plasma maximum radius
    real_t b = 1; //Placeholder: wall radius

    eqsys->SetUnknown(OptionConstants::UQTY_I_P, eqsys->GetScalarGrid());
    eqsys->SetUnknown(OptionConstants::UQTY_I_WALL, eqsys->GetScalarGrid());
    eqsys->SetUnknown(OptionConstants::UQTY_PSI_WALL, eqsys->GetScalarGrid());

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    const len_t id_psi_w = unknowns->GetUnknownID(OptionConstants::UQTY_PSI_WALL);
    const len_t id_psi_p = unknowns->GetUnknownID(OptionConstants::UQTY_POL_FLUX);
    const len_t id_I_w   = unknowns->GetUnknownID(OptionConstants::UQTY_I_WALL);
    const len_t id_I_p   = unknowns->GetUnknownID(OptionConstants::UQTY_I_P);
    const len_t id_j_tot = unknowns->GetUnknownID(OptionConstants::UQTY_J_TOT);

    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::Grid *fluidGrid  = eqsys->GetFluidGrid();
    FVM::Operator *eqn_pw1 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_pw2 = new FVM::Operator(scalarGrid);

//    FVM::TransientTerm *tt = new FVM::TransientTerm(scalarGrid,id_psi_w);
//    eqn_pw1->AddTerm(tt);
    eqn_pw1->AddTerm(new FVM::TransientTerm(scalarGrid,id_psi_w));
    eqn_pw2->AddTerm(new FVM::IdentityTerm(scalarGrid,R_W));
    
    eqsys->SetOperator(id_psi_w, id_psi_w, eqn_pw1, "dpsi_w/dt = R_w*I_w");
    eqsys->SetOperator(id_psi_w, id_I_w, eqn_pw2);

    FVM::Operator *eqn_Iw1 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Iw2 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Iw3 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Iw4 = new FVM::Operator(scalarGrid);

    eqn_Iw1->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqn_Iw2->AddTerm(new FVM::IdentityTerm(scalarGrid,-L_W));
    eqn_Iw3->AddTerm(new PoloidalFluxAtEdgeTerm(scalarGrid,fluidGrid,unknowns,id_psi_p));
    eqn_Iw4->AddTerm(new SOLMutualInductanceTerm(scalarGrid,a,b));

    eqsys->SetOperator(id_I_w, id_psi_w, eqn_Iw1, "psi_w = L_w*I_w + M_wp*I_p");
    eqsys->SetOperator(id_I_w, id_I_w,   eqn_Iw2);
    eqsys->SetOperator(id_I_w, id_psi_p, eqn_Iw3);
    eqsys->SetOperator(id_I_w, id_I_p,   eqn_Iw4);

    FVM::Operator *eqn_Ip1 = new FVM::Operator(scalarGrid);
    FVM::Operator *eqn_Ip2 = new FVM::Operator(scalarGrid);
    
    eqn_Ip1->AddTerm(new FVM::IdentityTerm(scalarGrid));
    eqn_Ip2->AddTerm(new TotalPlasmaCurrentFromJTot(scalarGrid,fluidGrid,unknowns,id_j_tot));

    eqsys->SetOperator(id_I_p, id_I_p,   eqn_Ip1, "Ip = integral(j_tot)");
    eqsys->SetOperator(id_I_p, id_j_tot, eqn_Ip2);

    eqsys->initializer->AddRule(
        id_psi_w,
        EqsysInitializer::INITRULE_EVAL_FUNCTION,
        [](FVM::UnknownQuantityHandler*, real_t *x){x[0] = 0;}
    );

    eqsys->initializer->AddRule(
        id_I_w,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_psi_w,
        id_I_p
    );

    eqsys->initializer->AddRule(
        id_I_p,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_j_tot
    );


}


