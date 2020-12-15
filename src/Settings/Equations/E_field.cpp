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
#include "FVM/Equation/DiagonalLinearTerm.hpp"
#include "FVM/Equation/LinearTransientTerm.hpp"
#include "DREAM/Equations/Fluid/HyperresistiveDiffusionTerm.hpp"

#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

/**
 * Implementation of a class which represents the Vloop term of the electric field diffusion equation.
 */
namespace DREAM {
    class VloopTerm : public FVM::DiagonalLinearTerm {
    public:
        VloopTerm(FVM::Grid* g) : FVM::DiagonalLinearTerm(g){}

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = 2*M_PI*grid->GetVpVol(ir)*sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}


/**
 * Implementation of a class which represents the dpsi/dt term of the electric field diffusion equation.
 */
namespace DREAM {
    class dPsiDtTerm : public FVM::LinearTransientTerm {
    public:
        dPsiDtTerm(FVM::Grid* g, const len_t unknownId) : FVM::LinearTransientTerm(g,unknownId){}

        virtual void SetWeights() override {
            len_t offset = 0;
            FVM::RadialGrid *rGrid = grid->GetRadialGrid();
            for (len_t ir = 0; ir < nr; ir++){
                real_t w = - rGrid->GetVpVol(ir)*rGrid->GetFSA_1OverR2(ir) * rGrid->GetBTorG(ir) / rGrid->GetBmin(ir);
                for(len_t i = 0; i < n1[ir]*n2[ir]; i++)
                    weights[offset + i] = w;
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
    

    // Type of boundary condition on the wall
    s->DefineSetting(MODULENAME "/bc/type", "Type of boundary condition to use on the wall for self-consistent E-field", (int_t)OptionConstants::UQTY_V_LOOP_WALL_EQN_SELFCONSISTENT);

    // Minor radius of the wall, defaults to radius of the plasma.
    s->DefineSetting(MODULENAME "/bc/wall_radius", "Minor radius of the inner wall", (real_t) -1);

    // Inverse wall time, defaults to 0 (infinitely conducting wall, 
    // which is equivalent to prescribing V_loop_wall to 0)
    s->DefineSetting(MODULENAME "/bc/inverse_wall_time", "Inverse wall time, representing the conductivity of the first wall", (real_t) 0.0);
    s->DefineSetting(MODULENAME "/bc/R0", "Major radius used to evaluate the external inductance for conductivity of the first wall", (real_t) 0.0);

    // TODO: Prescribed data (in time)
    DefineDataT(MODULENAME "/bc", s, "V_loop_wall");

    // Transport settings
    DefineOptions_Transport("eqsys/psi_p", s, false);

    
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
    FVM::Operator *eqn = new FVM::Operator(eqsys->GetFluidGrid());

    FVM::Interpolator1D *interp = LoadDataRT_intp(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s);
    FVM::PrescribedParameter *pp = new FVM::PrescribedParameter(eqsys->GetFluidGrid(), interp);
    eqn->AddTerm(pp);

    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, eqn, "Prescribed");
    // Initial value
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_E_FIELD,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );

    // Set boundary condition psi_wall = 0
    ConstructEquation_psi_wall_zero(eqsys,s);
}

/**
 * Construct the equation for a self-consistent electric field.
 */
void SimulationGenerator::ConstructEquation_E_field_selfconsistent(
    EquationSystem *eqsys, Settings* s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    // Set equations for self-consistent E field evolution
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

    // Add transient term -dpsi/dt
    Op1->AddTerm(new dPsiDtTerm(fluidGrid, eqsys->GetUnknownID(OptionConstants::UQTY_POL_FLUX)));
    // Add Vloop term
    Op2->AddTerm(new VloopTerm(fluidGrid));

    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op3, "eqsys/psi_p", fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI,
        eqsys->GetUnknownHandler(), s, false, false
    );

    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_POL_FLUX, Op1, "dpsi_p/dt = V_loop");
    eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_E_FIELD, Op2);
    
    if(hasTransport)
        eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, Op3, "dpsi_p/dt = V_loop + hyperresistivity(j_tot)");
    

    // for now: skip over the hyperresistive term
    bool setHyperresistivity = false;
    if(setHyperresistivity){
        // Add hyperresistivity term with placeholder values for Lambda and psi_t
        real_t *Lambda = new real_t[fluidGrid->GetNr()];
        real_t *psi_t = new real_t[fluidGrid->GetNr()];
        for(len_t ir=0; ir<fluidGrid->GetNr(); ir++){
            Lambda[ir] = 1;
            psi_t[ir]  = M_PI*fluidGrid->GetRadialGrid()->GetR(ir)*fluidGrid->GetRadialGrid()->GetR(ir)
                        * fluidGrid->GetRadialGrid()->GetBTorG(ir);
        }
        
        FVM::Operator *Op3 = new FVM::Operator(fluidGrid);
        Op3->AddTerm(new HyperresistiveDiffusionTerm(
        fluidGrid, Lambda, psi_t) 
        );
        eqsys->SetOperator(OptionConstants::UQTY_E_FIELD, OptionConstants::UQTY_J_TOT, Op3);

    } 

    /**
     * Load initial electric field profile.
     * If the input profile is not explicitly set, then 'SetInitialValue()' is
     * called with a null-pointer which results in E=0 at t=0
     */
    real_t *Efield_init = LoadDataR(MODULENAME, eqsys->GetFluidGrid()->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(OptionConstants::UQTY_E_FIELD, Efield_init);
    delete [] Efield_init;

    // Set equation for self-consistent boundary condition
    ConstructEquation_psi_wall_selfconsistent(eqsys,s);
}

