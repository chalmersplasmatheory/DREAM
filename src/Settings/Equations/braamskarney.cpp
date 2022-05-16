#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/EquationSystem.hpp"
#include "FVM/config.h"
#include "FVM/Equation/Operator.hpp"
#include "FVM/Equation/IdentityTerm.hpp"
#include "FVM/Equation/TransientTerm.hpp"
#include "FVM/Equation/MomentQuantity.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyPotentialIntegral.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyRHS.hpp"
#include "DREAM/Equations/Kinetic/BraamsKarneyLOperator.hpp"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>

using namespace DREAM;

// Integrands if neglecting xi-dependence.
static real_t IntegrandU0(real_t p, real_t pprime) {
    real_t a = sqrt((1 + p * p) * (1 + pprime * pprime));
	real_t b = p * pprime;
	
    if (b == 0)
        return -0.5 * 2.0 / sqrt(a * a - 1);
    
    real_t A = sqrt(a * a - 1);
    real_t B = -2 * a * b;
    real_t C = b * b;

    real_t low_lim = A - sqrt(A * A - B + C);
	real_t upp_lim = -A;
    
    if (A * A + B + C > 0) {
		upp_lim = sqrt(A * A + B + C) - A;
	}
	
	return -0.5 * (log((upp_lim + b) * (low_lim - b) / (upp_lim - b) / (low_lim + b))) / b;
}

static real_t IntegrandU1(real_t p, real_t pprime) {
	real_t a = sqrt((1 + p * p) * (1 + pprime * pprime));
	real_t b = p * pprime;
	
	if (b == 0)
		return -1.0/4.0 * 2 * sqrt(a * a - 1);
	
	auto indef = [] (real_t x) {return 1.0/2.0 * x * sqrt(x * x - 1) - 1.0/2.0 * log(sqrt(x * x - 1) + x);};
	real_t low = a - b;
	real_t high = a + b;
	if (low < 1)
		low = 1;
	
	return -1.0/4.0 * (indef(high) - indef(low)) / b;
}

static real_t IntegrandU2_arccosh_part(real_t p, real_t pprime) {
	real_t a = sqrt((1 + p * p) * (1 + pprime * pprime));
	real_t b = p * pprime;

    if (b == 0)
        return -1.0/16.0 * 2 * a * acosh(a);
    
    auto indef = [] (real_t x) { return 1.0/4.0 * (2 * x * x * acosh(x) - sqrt(x - 1) * sqrt(x + 1)*x - log(x + sqrt(x - 1) * sqrt(x + 1))); };
    real_t low = a - b;
    real_t high = a + b;
    if (low < 1)
        low = 1;
    
    return -1.0/16.0 * (indef(high) - indef(low)) / b;
}

static real_t IntegrandU2(real_t p, real_t pprime) {
	return IntegrandU2_arccosh_part(p, pprime) - IntegrandU1(p, pprime) / 4.0;
}

static real_t IntegrandPI0(real_t p, real_t pprime) {
    real_t a = sqrt((1 + p * p) * (1 + pprime * pprime));
	real_t b = p * pprime;
	
    if (b == 0)
        return -0.5 * 2.0 * a / sqrt(a * a - 1);
    
	auto indef = [] (real_t x) { return 1 / (x + 1) + 1 / (x - 1); };

	real_t indef_low = 0;
	if (a - b - 1 > 0) {
		real_t low_lim = sqrt(2 / (a - b - 1) + 1);
		indef_low = indef(low_lim);
	}

    real_t upp_lim = sqrt(2 / (a + b - 1) + 1);
	real_t indef_high = indef(upp_lim);
	
	return -0.5 * (indef_high - indef_low) / b;
}

static real_t IntegrandPI1(real_t p, real_t pprime) {
    real_t a = sqrt((1 + p * p) * (1 + pprime * pprime));
	real_t b = p * pprime;
	
    if (b == 0)
        return -1 / 4.0 * 2.0 * acosh(a);
    
    real_t low_lim = a - b;
    real_t upp_lim = a + b;

	if (low_lim < 1)
		low_lim = 1;
    
    auto indef = [] (real_t x) { return x * log(x + sqrt(x - 1) * sqrt(x + 1)) - sqrt(x - 1) * sqrt(x + 1); };
    
    return -1.0 / 4.0 * (indef(upp_lim) - indef(low_lim)) / b;
}

void InitializePotential(EquationSystem *eqsys, const char *unknown_name, const char *rhs_name, real_t a, real_t *f,
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

	auto *mg = hottailGrid->GetMomentumGrid(0);
	len_t np1 = mg->GetNp1();
	len_t np2 = mg->GetNp2();
	const real_t *p = mg->GetP1();
	const real_t *gamma = mg->GetGamma();
	const real_t *dp = mg->GetDp1();
	real_t *integrated = new real_t[np1];

    // FVM::Operator *eqn = new FVM::Operator(hottailGrid);
    // eqn->AddTerm(new FVM::TransientTerm(hottailGrid, id_pot));
    // eqsys->SetOperator(id_pot, id_pot, eqn);

	auto integrate = [np1, gamma, dp, p, f](real_t p_i, auto integrand) {
		real_t total = 0;
		for (len_t i = 0; i < np1; i++) {
			total += integrand(p_i, p[i]) * p[i] * p[i] / gamma[i] * f[i] * dp[i];
		}
		return total;
	};

	for (len_t i = 0; i < np1; i++) {
		switch (integralType) {
		case BraamsKarneyPotentialIntegral::Type::UPS0:
			integrated[i] = integrate(p[i], IntegrandU0);
			break;

		case BraamsKarneyPotentialIntegral::Type::UPS1:
			integrated[i] = integrate(p[i], IntegrandU1);
			break;

		case BraamsKarneyPotentialIntegral::Type::UPS2:
			integrated[i] = integrate(p[i], IntegrandU2);
			break;

		case BraamsKarneyPotentialIntegral::Type::PI0:
			integrated[i] = integrate(p[i], IntegrandPI0);
			break;

		case BraamsKarneyPotentialIntegral::Type::PI1:
			integrated[i] = integrate(p[i], IntegrandPI1);
			break;
		}
	}

	len_t nr = hottailGrid->GetNr();
    real_t *init = new real_t[hottailGrid->GetNCells()];

	len_t offset = 0;
	for (len_t ir = 0; ir < nr; ir++) {
		real_t *f = init + offset;
		for (len_t j = 0; j < np2; j++) {
			for (len_t i = 0; i < np1; i++) {
				f[j * np1 + i] = integrated[i];
			}
		}

		offset += np2 * np1;
	}

    eqsys->SetInitialValue(id_pot, init);

    // For a start set the initial value to 0.
	delete [] integrated;
	delete [] init;
}

void SimulationGenerator::ConstructEquation_BraamsKarney(EquationSystem *eqsys, Settings *s) {
    len_t nx[3];
    if (s->GetRealArray("eqsys/f_hot" "/init/x", 3, nx, false) != nullptr) {
		throw SettingsException("Braams-Karney can only be initialized to Maxwellian.");
	}

    FVM::Grid *hottailGrid = eqsys->GetHotTailGrid();

	real_t *n0 = LoadDataR("eqsys/f_hot", hottailGrid->GetRadialGrid(), s, "n0");
	real_t *T0 = LoadDataR("eqsys/f_hot", hottailGrid->GetRadialGrid(), s, "T0");

	// Momentum grid is assumed to be independant of radial coordinate.
	auto *mg = hottailGrid->GetMomentumGrid(0);
	len_t np1 = mg->GetNp1();
	const real_t *p = mg->GetP();
	real_t *f = new real_t[np1];

	for (len_t i = 0; i < np1; i++)
		f[i] = Constants::RelativisticMaxwellian(p[i], n0[0], T0[0]);

	InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_0_HOT, OptionConstants::UQTY_F_HOT, 0, f,
                        BraamsKarneyPotentialIntegral::Type::UPS0);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_1_HOT, OptionConstants::UQTY_POT_UPS_0_HOT, 2, f,
                        BraamsKarneyPotentialIntegral::Type::UPS1);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_UPS_2_HOT, OptionConstants::UQTY_POT_UPS_1_HOT, 2, f,
                        BraamsKarneyPotentialIntegral::Type::UPS2);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_PI_0_HOT, OptionConstants::UQTY_F_HOT, 1, f,
                        BraamsKarneyPotentialIntegral::Type::PI0);
    InitializePotential(eqsys, OptionConstants::UQTY_POT_PI_1_HOT, OptionConstants::UQTY_POT_PI_0_HOT, 1, f,
                        BraamsKarneyPotentialIntegral::Type::PI1);

	delete [] f;
	delete [] T0;
	delete [] n0;
}
