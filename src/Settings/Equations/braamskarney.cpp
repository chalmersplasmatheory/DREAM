#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "FVM/config.h"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/MomentQuantity.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyPotentialIntegral.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyRHS.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyLOperator.hpp"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>

using namespace DREAM;

void InitializePotential(EquationSystem *eqsys, const char *unknown_name, const char *rhs_name, real_t a,
                         BraamsKarneyPotentialIntegral::Type integralType) {
    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

    auto uh = eqsys->GetUnknownHandler();
    len_t id_pot = uh->GetUnknownID(unknown_name);
    len_t id_rhs = uh->GetUnknownID(rhs_name);
    len_t id_f_hot = uh->GetUnknownID(OptionConstants::UQTY_F_HOT);

    FVM::Operator *L_op = new FVM::Operator(hottailGrid);
    L_op->AddTerm(new BraamsKarneyLOperator(hottailGrid, a));

    auto rhsTerm = new BraamsKarneyRHS(hottailGrid);
    auto integralTerm = new BraamsKarneyPotentialIntegral(hottailGrid, integralType);

    if (id_f_hot == id_rhs) {
        FVM::Operator *rhs_integral_op = new FVM::Operator(hottailGrid);
        rhs_integral_op->AddTerm(rhsTerm);
        rhs_integral_op->AddTerm(integralTerm);

        eqsys->SetOperator(id_pot, id_rhs, rhs_integral_op);
    } else {
        FVM::Operator *rhs_op = new FVM::Operator(hottailGrid);
        rhs_op->AddTerm(rhsTerm);

        FVM::Operator *integral_op = new FVM::Operator(hottailGrid);
        integral_op->AddTerm(integralTerm);

        eqsys->SetOperator(id_pot, id_rhs, rhs_op);
        eqsys->SetOperator(id_pot, id_f_hot, integral_op);
    }

    eqsys->SetOperator(id_pot, id_pot, L_op);

    // For a start set the initial value to 0.
    eqsys->SetInitialValue(id_pot, nullptr);
}

void SimulationGenerator::ConstructEquation_BraamsKarney(EquationSystem *eqsys, Settings *) {
    InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_0_HOT, OptionConstants::UQTY_F_HOT, 0,
                        BraamsKarneyPotentialIntegral::Type::UPS0);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_1_HOT, OptionConstants::UQTY_POT_UPS_0_HOT, 2,
                        BraamsKarneyPotentialIntegral::Type::UPS1);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_2_HOT, OptionConstants::UQTY_POT_UPS_1_HOT, 2,
                        BraamsKarneyPotentialIntegral::Type::UPS2);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_PI_0_HOT, OptionConstants::UQTY_F_HOT, 1,
                        BraamsKarneyPotentialIntegral::Type::PI0);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_PI_1_HOT, OptionConstants::UQTY_POT_PI_0_HOT, 1,
                        BraamsKarneyPotentialIntegral::Type::PI1);
}
