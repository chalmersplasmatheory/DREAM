/**
 * Definition of equations relating to the electric field.
 *
 * Note that we define the electric field as
 *
 *      <E*B>
 *   -----------
 *   sqrt(<B^2>)
 *
 * in DREAM, where 'E' denotes the local electric field, 'B' the
 * local magnetic field and '<X>' denotes the flux-surface average
 * of a quantity 'X'.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Equation/WeightedIdentityTerm.hpp"
#include "FVM/Equation/WeightedTransientTerm.hpp"
#include "DREAM/Equations/PoloidalFlux/HyperresistiveDiffusionTerm.hpp"

#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

/**
 * Implementation of a class which represents the Vloop term of the electric field diffusion equation.
 */
namespace DREAM {
    class VloopTerm : public FVM::WeightedIdentityTerm {
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
    public:
        VloopTerm(FVM::Grid* g) : FVM::WeightedIdentityTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = 2*M_PI*sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}


/**
 * Implementation of a class which represents the dpsi/dt term of the electric field diffusion equation.
 */
namespace DREAM {
    class dPsiDtTerm : public FVM::WeightedTransientTerm {
    protected:
        virtual bool TermDependsOnUnknowns() override {return false;}
    public:
        dPsiDtTerm(FVM::Grid* g, const len_t unknownId) : FVM::WeightedTransientTerm(g,unknownId){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = - grid->GetRadialGrid()->GetFSA_1OverR2(ir) * grid->GetRadialGrid()->GetBTorG(ir) / grid->GetRadialGrid()->GetBmin(ir);
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}


#define MODULENAME "eqsys/E_field"


/**
 * Define options for the electric field module.
 */
void SimulationGenerator::DefineOptions_ElectricField(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electric field evolution", (int_t)OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED);

    // Prescribed data (in radius+time)
    DefineDataRT(MODULENAME, s, "data");

    // Prescribed initial profile (when evolving E self-consistently)
    DefineDataR(MODULENAME, s, "init");
    
}

/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_E_field(
    EquationSystem *eqsys, Settings *s
) {
    enum OptionConstants::uqty_E_field_eqn type = (enum OptionConstants::uqty_E_field_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_E_FIELD_EQN_PRESCRIBED:
            ConstructEquation_E_field_prescribed(eqsys, s);
            break;

        case OptionConstants::UQTY_E_FIELD_EQN_SELFCONSISTENT:
            ConstructEquation_E_field_selfconsistent(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_E_FIELD, type
            );
    }
}

/**
 * Construct the equation for a prescribed electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    // const real_t t0 = 0;
    FVM::Equation *eqn = new FVM::Equation(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn, "Prescribed");
    //eqsys->SetInitialValue(OptionConstants::UQTY_E_FIELD, interp->Eval(t0), t0);
}

/**
 * Construct the equation for a self-consistent electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_selfconsistent(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    // The self-consistent electric field requires an additional equation for the poloidal flux
    eqsys->SetUnknown(OptionConstants::UQTY_POL_FLUX, fluidGrid);
    ConstructEquation_psi_p(eqsys,s);

    // Set equations for self-consistent E field evolution
    FVM::Equation *eqn_E1 = new FVM::Equation(fluidGrid);
    FVM::Equation *eqn_E2 = new FVM::Equation(fluidGrid);

    // Add transient term -dpsi/dt
    eqn_E1->AddTerm(new dPsiDtTerm(fluidGrid, eqsys->GetUnknownID(OptionConstants::UQTY_POL_FLUX)));
    // Add Vloop term
    eqn_E2->AddTerm(new VloopTerm(fluidGrid));

    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_POL_FLUX, eqn_E1, "Poloidal flux resistive diffusion equation");
    eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn_E2);
    
    // for now: skip over the if statement
    bool settingHyperresistivity = false;
    if(settingHyperresistivity){
        // Add hyperresistivity term with placeholder values for Lambda and psi_t
        real_t *Lambda = new real_t[fluidGrid->GetNr()];
        real_t *psi_t = new real_t[fluidGrid->GetNr()];
        for(len_t ir=0; ir<fluidGrid->GetNr(); ir++){
            Lambda[ir] = 1;
            psi_t[ir]  = M_PI*fluidGrid->GetRadialGrid()->GetR(ir)*fluidGrid->GetRadialGrid()->GetR(ir)
                        * fluidGrid->GetRadialGrid()->GetBTorG(ir);
        }
        
        FVM::Equation *eqn_E3 = new FVM::Equation(fluidGrid);
        eqn_E3->AddTerm(new HyperresistiveDiffusionTerm(
        fluidGrid, Lambda, psi_t) 
        );
        eqsys->SetEquation(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, eqn_E3);


        real_t *Efield_init = LoadDataR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
        eqsys->SetInitialValue(OptionConstants::UQTY_E_FIELD, Efield_init);
        delete [] Efield_init;

    } 
}

