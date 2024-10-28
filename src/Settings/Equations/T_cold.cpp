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
#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"
#include "DREAM/Equations/Fluid/IonisationHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
#include "DREAM/Equations/Fluid/RadiatedPowerTerm.hpp"
#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIIonizLossTerm.hpp"
#include "DREAM/Equations/Fluid/ElectronHeatTerm.hpp"
#include "FVM/Equation/PrescribedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Fluid/ParallelHeatLossTerm.hpp"


using namespace DREAM;


#define MODULENAME "eqsys/T_cold"
#define MODULENAME_SPI "eqsys/spi"
#define MODULENAME_ION "eqsys/n_i"
#define MODULENAME_NRE "eqsys/n_re"


/**
 * Define options for the electron temperature module.
 */
void SimulationGenerator::DefineOptions_T_cold(Settings *s){
    s->DefineSetting(MODULENAME "/type", "Type of equation to use for determining the electron temperature evolution", (int_t)OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED);
    s->DefineSetting(MODULENAME "/recombination", "Whether to include recombination radiation (true) or ionization energy loss (false)", (bool)false);
    s->DefineSetting(MODULENAME "/parallel_losses", "Whether to include parallel losses (true) or not (false)", (bool)false);
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
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    enum OptionConstants::uqty_T_cold_eqn type = (enum OptionConstants::uqty_T_cold_eqn)s->GetInteger(MODULENAME "/type");

    switch (type) {
        case OptionConstants::UQTY_T_COLD_EQN_PRESCRIBED:
            ConstructEquation_T_cold_prescribed(eqsys, s);
            break;

        case OptionConstants::UQTY_T_COLD_SELF_CONSISTENT:
            ConstructEquation_T_cold_selfconsistent(eqsys, s, adas, nist, amjuel, oqty_terms);
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
    
    ConstructEquation_W_cold(eqsys, s);
}


/**
 * Construct the equation for a self-consistent temperature evolution.
 */
void SimulationGenerator::ConstructEquation_T_cold_selfconsistent(
    EquationSystem *eqsys, Settings *s, ADAS *adas, NIST *nist, AMJUEL *amjuel,
    struct OtherQuantityHandler::eqn_terms *oqty_terms
) {
    FVM::Grid *fluidGrid = eqsys->GetFluidGrid();

    IonHandler *ionHandler = eqsys->GetIonHandler();

    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    len_t id_T_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    len_t id_W_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_W_COLD);
    len_t id_n_cold  = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    len_t id_E_field = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);

    
    FVM::Operator *Op1 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op2 = new FVM::Operator(fluidGrid);
    FVM::Operator *Op3 = new FVM::Operator(fluidGrid);

    Op1->AddTerm(new FVM::TransientTerm(fluidGrid,id_W_cold) );

    // Check if parallel losses should be included
    bool lcfs_user_input_psi = (len_t)s->GetInteger(MODULENAME_NRE  "/lcfs_user_input_psi");
	real_t lcfs_psi_edge_t0 = s->GetReal(MODULENAME_NRE "/lcfs_psi_edge_t0");

    bool parallel_losses = s->GetBool(MODULENAME "/parallel_losses");
    if (parallel_losses) {
        Op1->AddTerm(new ParallelHeatLossTerm(fluidGrid,unknowns,ionHandler,-1,lcfs_user_input_psi, lcfs_psi_edge_t0)); // Add the term for parallel losses
    }


    oqty_terms->T_cold_ohmic = new OhmicHeatingTerm(fluidGrid,unknowns);
    Op2->AddTerm(oqty_terms->T_cold_ohmic);

    bool withRecombinationRadiation = s->GetBool(MODULENAME "/recombination");
    
    // Load opacity settings
    len_t nitypes;
    const int_t *iopacity_modes = s->GetIntegerArray(MODULENAME_ION "/opacity_modes", 1, &nitypes);
    enum OptionConstants::ion_opacity_mode *opacity_mode = new enum OptionConstants::ion_opacity_mode[nitypes];
    for (len_t i = 0; i < nitypes; i++)
        opacity_mode[i] = (enum OptionConstants::ion_opacity_mode)iopacity_modes[i];
        
    oqty_terms->T_cold_radiation = new RadiatedPowerTerm(
        fluidGrid,unknowns,ionHandler,adas,nist,amjuel,opacity_mode,withRecombinationRadiation
    );
    Op3->AddTerm(oqty_terms->T_cold_radiation);


    FVM::Operator *Op4 = new FVM::Operator(fluidGrid);
    // Add transport terms, if enabled
    bool hasTransport = ConstructTransportTerm(
        Op4, MODULENAME, fluidGrid,
        OptionConstants::MOMENTUMGRID_TYPE_PXI,
        eqsys, s, false, true,
        &oqty_terms->T_cold_advective_bc,&oqty_terms->T_cold_diffusive_bc
    );

    eqsys->SetOperator(id_T_cold, id_E_field,Op2);
    eqsys->SetOperator(id_T_cold, id_n_cold,Op3);
    std::string desc = "dWc/dt = j_ohm*E - sum_i n_cold*n_i*L_i";


    // SPI heat absorbtion
    OptionConstants::eqterm_spi_heat_absorbtion_mode spi_heat_absorbtion_mode = (enum OptionConstants::eqterm_spi_heat_absorbtion_mode)s->GetInteger(MODULENAME_SPI "/heatAbsorbtion");
    if(spi_heat_absorbtion_mode!=OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT){
        Op4->AddTerm(new SPIHeatAbsorbtionTerm(fluidGrid,eqsys->GetSPIHandler(),-1));
    }

    // SPI ionization losses (from neutral to equilibrium)
    
    // TODO: simplify the bool logic below if possible
    OptionConstants::eqterm_ionization_mode ionization_mode = 
        (enum OptionConstants::eqterm_ionization_mode)s->GetInteger(MODULENAME_ION "/ionization");
    bool includeKineticIonization = (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC) 
                                 || (ionization_mode == OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC);
    if(includeKineticIonization && !(eqsys->HasHotTailGrid()||eqsys->HasRunawayGrid()))
        throw SettingsException("Invalid ionization mode: cannot use kinetic ionization without a kinetic grid.");
    bool collfreqModeIsFull = (enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode")
        == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL;
    bool addFluidIonization = !(includeKineticIonization && eqsys->HasHotTailGrid() && collfreqModeIsFull);
    bool addFluidJacobian = (includeKineticIonization && eqsys->HasHotTailGrid() && (ionization_mode==OptionConstants::EQTERM_IONIZATION_MODE_KINETIC_APPROX_JAC));
    
    len_t nZSPInShard;
    OptionConstants::eqterm_spi_deposition_mode spi_deposition_mode = (enum OptionConstants::eqterm_spi_deposition_mode)s->GetInteger(MODULENAME_SPI "/deposition");
    OptionConstants::eqterm_spi_abl_ioniz_mode spi_abl_ioniz_mode = (enum OptionConstants::eqterm_spi_abl_ioniz_mode)s->GetInteger(MODULENAME_SPI "/abl_ioniz"); 
    const real_t *SPIMolarFraction  = s->GetRealArray(MODULENAME_ION "/SPIMolarFraction", 1, &nZSPInShard);
    
    // Add one ionization loss term for every species the pellet consists of
    // Note that, when accounting for the heat absorbed in the neutral cloud,
    // this energy is currently being redeposited where the material is deposited,
    // and we therefore need these terms even when heat absorption is included.
    if(spi_deposition_mode!=OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT && spi_abl_ioniz_mode!=OptionConstants::EQTERM_SPI_ABL_IONIZ_MODE_NEUTRAL){
        len_t offset=0;
        len_t nShard = eqsys->GetSPIHandler()->GetNShard();
        const len_t nZ = ionHandler->GetNZ();
        for(len_t iZ=0;iZ<nZ;iZ++){
            if(SPIMolarFraction[offset]>0){
                Op4->AddTerm(new IonSPIIonizLossTerm(fluidGrid, eqsys->GetIonHandler(), iZ, adas, eqsys->GetUnknownHandler(),
                    addFluidIonization, addFluidJacobian, eqsys->GetSPIHandler(), SPIMolarFraction,offset,1,nist,false, spi_abl_ioniz_mode));
                offset+=nShard;
            }else {
            	offset+=1;
            }
        }
    }

    if(hasTransport){
        oqty_terms->T_cold_transport = Op4->GetAdvectionDiffusion();
        desc += " + transport";
    }
    
    if (spi_deposition_mode!=OptionConstants::EQTERM_SPI_DEPOSITION_MODE_NEGLECT || 
        spi_heat_absorbtion_mode!=OptionConstants::EQTERM_SPI_HEAT_ABSORBTION_MODE_NEGLECT ||
        hasTransport)
        eqsys->SetOperator(id_T_cold, id_T_cold,Op4); 

    /**
     * If hot-tail grid is enabled, add collisional  
     * energy transfer from hot-tail to T_cold. 
     * NOTE: We temporarily disable the collisional energy transfer 
     *       in collfreq_mode FULL because we lack the correction 
     *       term for when electrons transfer from the cold to the 
     *       hot region. This should be corrected for at some point!
     */
    bool collfreqModeFull = ((enum OptionConstants::collqty_collfreq_mode)s->GetInteger("collisions/collfreq_mode") == OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_FULL);
    if( eqsys->HasHotTailGrid()&& !collfreqModeFull ){
        len_t id_f_hot = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);

        FVM::MomentQuantity::pThresholdMode pMode = 
            (FVM::MomentQuantity::pThresholdMode)s->GetInteger("eqsys/f_hot/pThresholdMode");
        real_t pThreshold = 0.0;
        if(collfreqModeFull){
            // With collfreq_mode FULL, only add contribution from hot electrons
            // defined as those with momentum above the defined threshold. 
            pThreshold = (real_t)s->GetReal("eqsys/f_hot/pThreshold");
        }
        oqty_terms->T_cold_fhot_coll = new CollisionalEnergyTransferKineticTerm(
            fluidGrid,eqsys->GetHotTailGrid(),
            id_T_cold, id_f_hot,eqsys->GetHotTailCollisionHandler(), eqsys->GetUnknownHandler(),
            eqsys->GetHotTailGridType(), -1.0,
            pThreshold, pMode
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( oqty_terms->T_cold_fhot_coll );
        eqsys->SetOperator(id_T_cold, id_f_hot, Op5);
        desc += " + int(W*nu_E*f_hot)";
    }
    // If runaway grid and not FULL collfreqmode, add collisional  
    // energy transfer from runaways to T_cold. 
    if( eqsys->HasRunawayGrid() ){
        len_t id_f_re = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);

        oqty_terms->T_cold_fre_coll = new CollisionalEnergyTransferKineticTerm(
            fluidGrid,eqsys->GetRunawayGrid(),
            id_T_cold, id_f_re,eqsys->GetRunawayCollisionHandler(),eqsys->GetUnknownHandler(),
            eqsys->GetRunawayGridType(), -1.0
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( oqty_terms->T_cold_fre_coll );
        eqsys->SetOperator(id_T_cold, id_f_re, Op5);
        desc += " + int(W*nu_E*f_re)";
    } else {
        len_t id_n_re = unknowns->GetUnknownID(OptionConstants::UQTY_N_RE);
        oqty_terms->T_cold_nre_coll = new CollisionalEnergyTransferREFluidTerm(
            fluidGrid, eqsys->GetUnknownHandler(), eqsys->GetREFluid()->GetLnLambda(), -1.0
        );
        FVM::Operator *Op5 = new FVM::Operator(fluidGrid);
        Op5->AddTerm( oqty_terms->T_cold_nre_coll );
        eqsys->SetOperator(id_T_cold, id_n_re, Op5);

        desc += " + e*c*Ec*n_re";
    }
    
    // ADD COLLISIONAL ENERGY TRANSFER WITH ION SPECIES
    OptionConstants::uqty_T_i_eqn Ti_type = 
        (OptionConstants::uqty_T_i_eqn)s->GetInteger("eqsys/n_i/typeTi");
    if(Ti_type == OptionConstants::UQTY_T_I_INCLUDE) {
        CoulombLogarithm *lnLambda = eqsys->GetREFluid()->GetLnLambda();
        const len_t nZ = ionHandler->GetNZ();
        const len_t id_Wi = eqsys->GetUnknownID(OptionConstants::UQTY_WI_ENER);
        oqty_terms->T_cold_ion_coll = new FVM::Operator(fluidGrid);
        for(len_t iz=0; iz<nZ; iz++){
            oqty_terms->T_cold_ion_coll->AddTerm(
                new MaxwellianCollisionalEnergyTransferTerm(
                    fluidGrid,
                    0, false,
                    iz, true,
                    unknowns, lnLambda, ionHandler, -1.0)
            );
        }
        eqsys->SetOperator(id_T_cold, id_Wi, oqty_terms->T_cold_ion_coll);
        desc += " + sum_i Q_ei";
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
    Op2->AddTerm(new ElectronHeatTerm(fluidGrid,id_n_cold,eqsys->GetUnknownHandler()) );

    eqsys->SetOperator(id_W_cold, id_W_cold, Op1, "W_cold = (3/2)*n_cold*T_cold");
    eqsys->SetOperator(id_W_cold, id_T_cold, Op2);    

    eqsys->initializer->AddRule(
        id_W_cold,
        EqsysInitializer::INITRULE_EVAL_EQUATION,
        nullptr,
        id_T_cold,
        id_n_cold
    );


}
