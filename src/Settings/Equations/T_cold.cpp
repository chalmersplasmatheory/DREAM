/**
 * Definition of equations relating to the cold electron temperature.
 */

#include "DREAM/EquationSystem.hpp"
#include "DREAM/NIST.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/IonisationHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/BindingEnergyTerm.hpp"
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


#define MODULENAME "eqsys/T_cold"


/**
 * Define options for the electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_cold(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);

    // Prescribed data (in radius+time)
    DefineDataRT(MODULENAME, s, "data");

    // Prescribed initial profile (when evolving T self-consistently)
    DefineDataR(MODULENAME, s, "init");
    
    // Transport settings
    DefineOptions_Transport(MODULENAME, s, false);
}


/**
 * Construct the equation for the electric field.
 */
void SimulationGenerator::ConstructEquation_T_cold(
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    enum OptionConstants::uqty_T_cold_eqn type = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        case OptionConstants::UQTY_T_COLD_SELF_CONSISTENT:
            ConstructEquation_T_cold_selfconsistent(eqsys, s, adas, nist, oqty_terms);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_T_COLD, type
            );
    }
}

/**
 * Construct the equation for a prescribed temperature.
 */
void SimulationGenerator::ConstructEquation_T_cold_prescribed(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    FVM::Operator *eqn = new FVM::Operator(fluidGrid);

    FVM::Interpolator1D *interp = LoadDataRT_intp(MODULENAME,fluidGrid->GetRadialGrid(), s);
    eqn->AddTerm(new FVM::PrescribedParameter(fluidGrid, interp));

    eqsys->SetOperator(OptionConstants::UQTY_T_COLD, OptionConstants::UQTY_T_COLD, eqn, "Prescribed");

    // Initialization
    eqsys->initializer->AddRule(
        OptionConstants::UQTY_T_COLD,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );
}


/**
 * Construct the equation for a self-consistent temperature evolution.
 */
void SimulationGenerator::ConstructEquation_T_cold_selfconsistent(
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    /**
     * The self-consistent temperature evolution uses an equation
     * for the total cold electron energy W_c (potential + heat) 
     */
    eqsys->SetUnknown(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD_DESC, fluidGrid);
    
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    len_t id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_W_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    len_t id_n_hot   = unknowns->GetUnknownID(OptionConstants::UQTY_N_HOT);
    len_t id_E_field = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::TransientTerm(fluidGrid,id_W_cold) );
    Op2->AddTerm(new OhmicHeatingTerm(fluidGrid,unknowns));
    oqty_terms->T_cold_radterm = new RadiatedPowerTerm(fluidGrid,unknowns,eqsys->GetIonHandler(),adas);
    Op3->AddTerm(oqty_terms->T_cold_radterm);


    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op4, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, s, false
    );

    std::string desc = "dWc/dt = j_ohm*E - sum_i n_cold*n_i*L_i";
    eqsys->SetOperator(id_T_cold, id_E_field,Op2);
    eqsys->SetOperator(id_T_cold, id_n_cold,Op3);

    if(hasTransport){
        eqsys->SetOperator(id_T_cold, id_T_cold,Op4);
        desc += " + transport";
    }
    // If hot-tail grid is enabled, add collisional  
    // energy transfer from hot-tail to T_cold. 
    if( eqsys->HasHotTailGrid() ){
        len_t id_f_hot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);

        FVM::MomentQuantity::pThresholdMode pMode = (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
        real_t pThreshold = 0.0;
        enum OptionConstants::collqty_collfreq_mode collfreq_mode =
            (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
        if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
            // With collfreq_mode FULL, only add contribution from hot electrons
            // defined as those with momentum above the defined threshold. 
            pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        }

        FVM::Operator *Op4 = new FVM::Operator(fluidGrid);


        Op4->AddTerm( new CollisionalEnergyTransferKineticTerm(fluidGrid,eqsys->GetHotTailGrid(),
            id_T_cold, id_f_hot,eqsys->GetHotTailCollisionHandler(), eqsys->GetUnknownHandler(), -1.0,
            pThreshold, pMode));
        eqsys->SetOperator(id_T_cold, id_f_hot, Op4);

        // TODO: IonisationHeatingTerm here is the old approximate hot-electron ionization correction. Should be replaced.
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( new IonisationHeatingTerm(fluidGrid, unknowns, eqsys->GetIonHandler(), adas, nist) );
        eqsys->SetOperator(id_T_cold, id_n_hot, Op5);
        desc += " - int(nu_E*f_hot)";
    }
    // If runaway grid and not FULL collfreqmode, add collisional  
    // energy transfer from runaways to T_cold. 
    if( eqsys->HasRunawayGrid() ){
        len_t id_f_re = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

        FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
        Op4->AddTerm( new CollisionalEnergyTransferKineticTerm(fluidGrid,eqsys->GetRunawayGrid(),
            id_T_cold, id_f_re,eqsys->GetRunawayCollisionHandler(),eqsys->GetUnknownHandler()));
        eqsys->SetOperator(id_T_cold, id_f_re, Op4);
        desc += " - int(nu_E*f_re)";
    }
    
    eqsys->SetOperator(id_T_cold, id_W_cold,Op1,desc);

    /**
     * Load initial electron temperature profile.
     * If the input profile is not explicitly set, then 'SetInitialValue()' is
     * called with a null-pointer which results in T=0 at t=0
     */
    real_t *Tcold_init = LoadDataR(MODULENAME, fluidGrid->GetRadialGrid(), s, "init");
    eqsys->SetInitialValue(id_T_cold, Tcold_init);
    delete [] Tcold_init;


    ConstructEquation_W_cold(eqsys, s, nist);
}


/**
 * Implementation of an equation term which represents the total
 * heat of the cold electrons: W_h = (3/2) * n_cold * T_cold
 */
namespace DREAM {
    class ElectronHeatTerm : public FVM::DiagonalQuadraticTerm {
    public:
        ElectronHeatTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u) 
            : FVM::DiagonalQuadraticTerm(g,u->GetUnknownID(OptionConstants::UQTY_N_COLD),u){}

        virtual void SetWeights() override {
            for(len_t i = 0; i<grid->GetNCells(); i++)
                weights[i] = 1.5 * Constants::ec;
        }
    };
}


/**
 * Construct the equation for electron energy content:
 *    W_cold = 3n_cold*T_cold/2 + W_binding,
 * where W_binding is the total binding energy of all
 * ions (i.e. the minimum energy required to fully ionise
 * the entire plasma). 
*/
void SimulationGenerator::ConstructEquation_W_cold(
    EquationSystem *eqsys, Settings* /*s*/, NIST* nist
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);

    len_t id_W_cold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_T_cold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_n_cold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    len_t id_n_i = eqsys->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0) );
    Op2->AddTerm(new ElectronHeatTerm(fluidGrid,eqsys->GetUnknownHandler()) );
    Op3->AddTerm(new BindingEnergyTerm(fluidGrid, eqsys->GetIonHandler(), nist));

    eqsys->SetOperator(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD, Op1, "W_c = 3nT/2 + W_bind");
    eqsys->SetOperator(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_T_COLD, Op2);    
    eqsys->SetOperator(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_ION_SPECIES, Op3);

    eqsys->initializer->AddRule(
        id_W_cold,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_T_cold,
        id_n_i,
        id_n_cold
    );


}
