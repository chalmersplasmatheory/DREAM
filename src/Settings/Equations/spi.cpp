#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/Fluid/SPIAblationTerm.hpp"
#include "DREAM/Equations/Fluid/IonSPIDepositionTerm.hpp"
#include "DREAM/Equations/Fluid/SPIHeatAbsorbtionTerm.hpp"
#include "DREAM/Equations/Fluid/SPITransientTerm.hpp"
#include "DREAM/Equations/Scalar/ConstantSPIVelocityTerm.hpp"
#include "DREAM/Equations/Scalar/ConstantSPIVelocityPositionTerm.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/ConstantParameter.hpp"
#include "FVM/Grid/Grid.hpp"


using namespace DREAM;

#define MODULENAME "eqsys/spi"

void SimulationGenerator::ConstructEquation_SPI(
    EquationSystem *eqsys, Settings *s
){
    enum OptionConstants::eqterm_spi_velocity_mode spi_velocity_mode = (enum OptionConstants::eqterm_spi_velocity_mode)s->GetInteger(MODULENAME "/velocity");
    enum OptionConstants::eqterm_spi_ablation_mode spi_ablation_mode = (enum OptionConstants::eqterm_spi_ablation_mode)s->GetInteger(MODULENAME "/ablation");

    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::UnknownQuantityHandler *unknowns = eqsys->GetUnknownHandler();

    // Get data for shard radii
    len_t nShard;
    const real_t *rp_init = s->GetRealArray(MODULENAME "/init/rp", 1, &nShard);

    /*eqsys->SetUnknown(OptionConstants::UQTY_R_P,scalarGrid,nShard);
    eqsys->SetUnknown(OptionConstants::UQTY_X_P,scalarGrid,3*nShard);
    eqsys->SetUnknown(OptionConstants::UQTY_V_P,scalarGrid,3*nShard);*/

    len_t id_rp=unknowns->GetUnknownID(OptionConstants::UQTY_R_P);
    //len_t id_vp=unknowns->GetUnknownID(OptionConstants::UQTY_V_P);
    //len_t id_xp=unknowns->GetUnknownID(OptionConstants::UQTY_X_P);

    // Ablation terms
    FVM::Operator *Op_rp = new FVM::Operator(scalarGrid);
    Op_rp->AddTerm(new SPITransientTerm(scalarGrid,id_rp,nShard) );

    if(spi_ablation_mode!=OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT)
        Op_rp->AddTerm(new SPIAblationTerm(scalarGrid, unknowns, eqsys->GetSPIHandler(), -1));

    eqsys->SetOperator(OptionConstants::UQTY_R_P, OptionConstants::UQTY_R_P, Op_rp);

    // Initialize shard radii
    eqsys->SetInitialValue(id_rp, rp_init);

    // Shard velocity and position terms
    switch (spi_velocity_mode) {
        case OptionConstants::EQTERM_SPI_VELOCITY_MODE_PRESCRIBED:
            ConstructEquation_v_p_prescribed_constant(eqsys, s);
            ConstructEquation_x_p_prescribed_constant_velocity(eqsys, s);
            break;

        default:
            throw SettingsException(
                "Unrecognized equation type for '%s': %d.",
                OptionConstants::UQTY_V_P, spi_velocity_mode
            );
    }
}

void SimulationGenerator::ConstructEquation_v_p_prescribed_constant(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::UnknownQuantityHandler* u = eqsys->GetUnknownHandler();
    FVM::Operator *eqn = new FVM::Operator(scalarGrid);

    // Initialize shard velocity
    len_t nShard;
    const real_t *vp_init;
    vp_init=s->GetRealArray(MODULENAME "/init/vp", 1, &nShard);

    eqn->AddTerm(new ConstantSPIVelocityTerm(scalarGrid, u, vp_init));

    eqsys->SetOperator(OptionConstants::UQTY_V_P, OptionConstants::UQTY_V_P, eqn, "Prescribed");

    // Initialization
    eqsys->SetInitialValue(OptionConstants::UQTY_V_P, vp_init);
    /*eqsys->initializer->AddRule(
        OptionConstants::UQTY_V_P,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );*/
}

void SimulationGenerator::ConstructEquation_x_p_prescribed_constant_velocity(
    EquationSystem *eqsys, Settings *s
) {
    FVM::Grid *scalarGrid = eqsys->GetScalarGrid();
    FVM::UnknownQuantityHandler* u = eqsys->GetUnknownHandler();
    FVM::Operator *eqn = new FVM::Operator(scalarGrid);

    // Initialize shard velocity
    len_t nShard;
    const real_t *xp_init;
    xp_init=s->GetRealArray(MODULENAME "/init/xp", 1, &nShard);

    eqn->AddTerm(new ConstantSPIVelocityPositionTerm(scalarGrid, u, xp_init));

    eqsys->SetOperator(OptionConstants::UQTY_X_P, OptionConstants::UQTY_X_P, eqn, "Prescribed");

    // Initialization
    eqsys->SetInitialValue(OptionConstants::UQTY_X_P, xp_init);
   /* eqsys->initializer->AddRule(
        OptionConstants::UQTY_X_P,
        EqsysInitializer::INITRULE_EVAL_EQUATION
    );*/
}
