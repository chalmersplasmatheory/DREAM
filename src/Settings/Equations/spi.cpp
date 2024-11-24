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

    // Get data for shard content
    len_t nShard;
    const real_t *rp_init = s->GetRealArray(MODULENAME "/init/rp", 1, &nShard);
    
    // Due to numerical advantages, we follow the shard radii using the variable Yp=rp^(5/3) in the code,
    // so calculate the initial values of this variable
    real_t *Yp_init = new real_t[nShard];
    for(len_t i=0;i<nShard;i++)
        Yp_init[i] = pow(rp_init[i],5.0/3.0);

    len_t id_Yp=unknowns->GetUnknownID(OptionConstants::UQTY_Y_P);

    // Ablation terms
    FVM::Operator *Op_Yp = new FVM::Operator(scalarGrid);
    Op_Yp->AddTerm(new SPITransientTerm(scalarGrid,id_Yp,nShard) );

    if(spi_ablation_mode!=OptionConstants::EQTERM_SPI_ABLATION_MODE_NEGLECT)
        Op_Yp->AddTerm(new SPIAblationTerm(scalarGrid, unknowns, eqsys->GetSPIHandler(), -1));

    eqsys->SetOperator(OptionConstants::UQTY_Y_P, OptionConstants::UQTY_Y_P, Op_Yp);

    // Initialize shard radii-variable
    eqsys->SetInitialValue(id_Yp, Yp_init);
	delete [] Yp_init;

    // Shard velocity and position terms
    switch (spi_velocity_mode) {
        case OptionConstants::EQTERM_SPI_VELOCITY_MODE_PRESCRIBED:
            ConstructEquation_v_p_prescribed_constant(eqsys, s);
            ConstructEquation_x_p_prescribed_constant_velocity(eqsys, s);
            break;
        case OptionConstants::EQTERM_SPI_VELOCITY_MODE_NONE:
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
    
    const real_t *t_delay;
    t_delay=s->GetRealArray(MODULENAME "/init/t_delay", 1, &nShard);
    
    eqn->AddTerm(new ConstantSPIVelocityPositionTerm(scalarGrid, u, xp_init, t_delay));

    eqsys->SetOperator(OptionConstants::UQTY_X_P, OptionConstants::UQTY_X_P, eqn, "Prescribed");

    // Initialization
    eqsys->SetInitialValue(OptionConstants::UQTY_X_P, xp_init);
}
