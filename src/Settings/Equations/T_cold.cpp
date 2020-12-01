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
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIIonizLossTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;


#define MODULENAME "eqsys/T_cold"
#define MODULENAME_SPI "eqsys/spi"
#define MODULENAME_ION "eqsys/n_i"


/**
 * Define options for the electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_cold(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);
    s->DefineSetting(MODULENAME "/recombination", "Whether to include recombination radiation (true) or ionization energy loss (false)", (bool)true);

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
    
    //eqsys->SetUnknown(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD_DESC, fluidGrid);
    ConstructEquation_W_cold(eqsys, s);
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
     * for the total cold electron heat energy W_c
     */
    //eqsys->SetUnknown(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD_DESC, fluidGrid);
    
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();
    len_t id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_W_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    len_t id_E_field = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::TransientTerm(fluidGrid,id_W_cold) );
    oqty_terms->T_cold_ohmic = new OhmicHeatingTerm(fluidGrid,unknowns);
    Op2->AddTerm(oqty_terms->T_cold_ohmic);

    bool withRecombinationRadiation = s->GetBool(MODULENAME "/recombination");
    oqty_terms->T_cold_radiation = new RadiatedPowerTerm(
        fluidGrid,unknowns,eqsys->GetIonHandler(),adas,nist,withRecombinationRadiation
    );
    Op3->AddTerm(oqty_terms->T_cold_radiation);


    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op4, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI,
        unknowns, s, false, true,
        &oqty_terms->T_cold_advective_bc,&oqty_terms->T_cold_diffusive_bc
    );

    eqsys->SetOperator(id_T_cold, id_E_field,Op2);
    eqsys->SetOperator(id_T_cold, id_n_cold,Op3);
    std::string desc = "dWc/dt = j_ohm*E - sum_i n_cold*n_i*L_i";


    // SPI heat absorbtion
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode = (enum OptionConstants::eqterm_spi_heat_absorbtion_mode)s->GetInteger(MODULENAME_SPI "/heatAbsorbtion");
    if(spi_heat_absorbtion_mode!=OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT){
        FVM::Operator *OpSPIHeat = new FVM::Operator(fluidGrid);
        OpSPIHeat->AddTerm(new SPIHeatAbsorbtionTerm(fluidGrid,eqsys->GetSPIHandler(),-1));
        eqsys->SetOperator(id_T_cold, id_T_cold, OpSPIHeat);
    }

    // SPI ionization losses
    len_t nZ;
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode = (enum OptionConstants::eqterm_spi_deposition_mode)s->GetInteger(MODULENAME_SPI "/deposition"); 
    const real_t *SPIMolarFraction  = s->GetRealArray(MODULENAME_ION "/SPIMolarFraction", 1, &nZ);
    if(spi_deposition_mode!=OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT){
        FVM::Operator *OpIonLoss = new FVM::Operator(fluidGrid);
        for(len_t iZ=0;iZ<nZ;iZ++){
            if(SPIMolarFraction[iZ]>0){
                OpIonLoss->AddTerm(new IonSPIIonizLossTerm(fluidGrid, eqsys->GetIonHandler(), iZ, eqsys->GetSPIHandler(), SPIMolarFraction[iZ],1,nist));
            }
        }
        eqsys->SetOperator(id_T_cold, id_T_cold, OpIonLoss);
    }

    if(hasTransport){
        oqty_terms->T_cold_transport = Op4->GetAdvectionDiffusion();
        eqsys->SetOperator(id_T_cold, id_T_cold,Op4);
        desc += " + transport";
    }
    // If hot-tail grid is enabled, add collisional  

    // energy transfer from hot-tail to T_cold. 
    if( eqsys->HasHotTailGrid() ){
        len_t id_f_hot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);

        FVM::MomentQuantity::pThresholdMode pMode = 
            (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
        real_t pThreshold = 0.0;
        enum OptionConstants::collqty_collfreq_mode collfreq_mode =
            (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode");
        if(collfreq_mode == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL){
            // With collfreq_mode FULL, only add contribution from hot electrons
            // defined as those with momentum above the defined threshold. 
            pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        }
        oqty_terms->T_cold_fhot_coll = new CollisionalEnergyTransferKineticTerm(
            fluidGrid,eqsys->GetHotTailGrid(),
            id_T_cold, id_f_hot,eqsys->GetHotTailCollisionHandler(), eqsys->GetUnknownHandler(), -1.0,
            pThreshold, pMode
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( oqty_terms->T_cold_fhot_coll );
        eqsys->SetOperator(id_T_cold, id_f_hot, Op5);
        desc += " + int(nu_E*f_hot)";
    }
    // If runaway grid and not FULL collfreqmode, add collisional  
    // energy transfer from runaways to T_cold. 
    if( eqsys->HasRunawayGrid() ){
        len_t id_f_re = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

        oqty_terms->T_cold_fre_coll = new CollisionalEnergyTransferKineticTerm(
            fluidGrid,eqsys->GetRunawayGrid(),
            id_T_cold, id_f_re,eqsys->GetRunawayCollisionHandler(),eqsys->GetUnknownHandler()
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( oqty_terms->T_cold_fre_coll );
        eqsys->SetOperator(id_T_cold, id_f_re, Op5);
        desc += " + int(nu_E*f_re)";
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


    ConstructEquation_W_cold(eqsys, s);
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
 *    W_cold = 3n_cold*T_cold/2
*/
void SimulationGenerator::ConstructEquation_W_cold(
    EquationSystem *eqsys, Settings* /*s*/
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();
    
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);

    len_t id_W_cold = eqsys->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_T_cold = eqsys->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_n_cold = eqsys->GetUnknownID(OptionConstants::UQTY_N_COLD);
    
    Op1->AddTerm(new FVM::IdentityTerm(fluidGrid,-1.0) );
    Op2->AddTerm(new ElectronHeatTerm(fluidGrid,eqsys->GetUnknownHandler()) );

    eqsys->SetOperator(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_W_COLD, Op1, "W_c = 3nT/2");
    eqsys->SetOperator(OptionConstants::UQTY_W_COLD, OptionConstants::UQTY_T_COLD, Op2);    

    eqsys->initializer->AddRule(
        id_W_cold,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_T_cold,
        id_n_cold
    );


}
