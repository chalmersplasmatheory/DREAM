/**
 * Equation system initializer.
 */

#include "DREAM/EqsysInitializer.hpp"
#include "DREAM/IO.hpp"


using namespace DREAM;
using namespace std;


/**
 * Constructor.
 */
EqsysInitializer::EqsysInitializer(
    FVM::UnknownQuantityHandler *unknowns,
    vector<UnknownQuantityEquation*> *unknown_equations
) : unknowns(unknowns), unknown_equations(unknown_equations) { }

/**
 * Destructor.
 */
EqsysInitializer::~EqsysInitializer() { }

/**
 * Add a new initialization rule.
 *
 * uqtyId: ID of the unknown quantity to which the rule applies.
 * type:   Rule type (see 'enum initrule_t').
 */
void EqsysInitializer::AddRule(const string& uqtyName, const enum initrule_t type, initfunc_t fnc) {
    const len_t uqtyId = this->unknowns->GetUnknownID(uqtyName);
    this->AddRule(uqtyId, type, fnc);
}
void EqsysInitializer::AddRule(const len_t uqtyId, const enum initrule_t type, initfunc_t fnc) {
    if (this->HasRuleFor(uqtyId))
        throw EqsysInitializerException(
            "An initialization rule for '%s' has already been added.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    if (type == INITRULE_STEADY_STATE_SOLVE)
        throw NotImplementedException(
            "%s: No support for solving for initial steady state implemented yet.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    this->rules[uqtyId] = new struct initrule(
        uqtyId, type, vector<len_t>(0), fnc
    );
}

/**
 * Execute all initialization rules and initialize the unknown
 * quantities of the equation system.
 *
 * t0: Time for which the system should be initialized.
 */
void EqsysInitializer::Execute(const real_t t0) {
    // List specifying order of execution
    vector<len_t> order;

    // Build the 'order' list
    for (auto rule : this->rules) {
        // Check if the rule has already been executed
        if (rule.second->executed)
            continue;

        vector<len_t> tmp = ConstructExecutionOrder(rule.second);
        order.insert(order.end(), tmp.begin(), tmp.end());
    }

    // Execute the rules in the calculated order
    for (len_t uqtyId : order) {
        struct initrule *rule = this->rules[uqtyId];

        switch (rule->type) {
            case INITRULE_EVAL_EQUATION:
                this->EvaluateEquation(t0, uqtyId);
                break;

            case INITRULE_EVAL_FUNCTION:
                this->EvaluateFunction(t0, uqtyId);
                break;

            case INITRULE_STEADY_STATE_SOLVE:
                throw NotImplementedException("The 'STEADY_STATE_SOLVE' initialization rule has not been implemented yet.");
                break;

            default:
                throw EqsysInitializerException(
                    "Unrecognized initialization rule type: %d.",
                    rule->type
                );
        }
    }

    // Verfiy that all unknown quantities have now been initialized
    this->VerifyAllInitialized();
}

/**
 * Calculates the order in which quantities on which the specified
 * rule depends should be initialized. The ID of the specified rule
 * is added last to the execution order.
 *
 * rule: Rule to resolve dependency execution order for.
 */
vector<len_t> EqsysInitializer::ConstructExecutionOrder(struct initrule *rule) {
    vector<len_t> order;

    // Traverse dependencies and ensure
    for (len_t dep : rule->dependencies) {
        if (this->unknowns->HasInitialValue(dep)) {
            continue;
        } else if (this->HasRuleFor(dep)) {
            struct initrule *deprule = this->rules[dep];

            // If the rule has already been executed, we can
            // safely proceed to the next dependency
            if (deprule->executed)
                continue;
            // Seeing this rule a second time means we have
            // discovered a circular dependency...
            else if (deprule->marked) {
                throw EqsysInitializerException(
                    "Unable to resolve circular dependency: %s -> %s.",
                    this->unknowns->GetUnknown(rule->uqtyId)->GetName().c_str(),
                    this->unknowns->GetUnknown(deprule->uqtyId)->GetName().c_str()
                );
            } else {
                // This rule has no un-initialized dependencies so we add
                // it to the list of rules.
                deprule->marked = true;
                vector<len_t> tmp = ConstructExecutionOrder(deprule);
                order.insert(order.end(), tmp.begin(), tmp.end());
            }
        } else
            throw EqsysInitializerException(
                "Unable to resolve initialization dependencies. No rule to initialize '%s'.",
                this->unknowns->GetUnknown(dep)->GetName()
            );
    }

    rule->executed = true;
    order.push_back(rule->uqtyId);
    return order;
}

/**
 * Check an initialization rule for the specified unknown
 * quantity exists.
 */
bool EqsysInitializer::HasRuleFor(const len_t uqtyId) const {
    return (this->rules.find(uqtyId) != this->rules.end());
}

/**
 * Verify that all unknown quantities in the UnknownQuantityHandler
 * have been successfully initialized.
 */
void EqsysInitializer::VerifyAllInitialized() const {
    bool initialized = true;
    for (auto it : this->unknowns->GetUnknowns()) {
        if (!it->HasInitialValue()) {
            initialized = false;
            DREAM::IO::PrintError(
                "%s: No initial value defined for unknown quantity.",
                it->GetName().c_str()
            );
        }
    }

    if (!initialized)
        throw EqsysInitializerException(
            "Some quantities have no defined initial values."
        );
}


/***********************
 * EVALUATION ROUTINES *
 ***********************/
/**
 * Initialize the specified unknown quantity by evaluating
 * its equation.
 *
 * t0:     Time for which the system should be initialized.
 * uqtyId: ID of unknown quantity to evaluate equation for.
 */
void EqsysInitializer::EvaluateEquation(const real_t t0, const len_t uqtyId) {
    UnknownQuantityEquation *eqn = this->unknown_equations->at(uqtyId);

    if (!eqn->IsEvaluable())
        throw EqsysInitializerException(
            "Unable to initialize '%s': equation is not evaluable.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );
    
    // Construct vector for temporary storage of quantity value
    const len_t N = eqn->NumberOfElements();
    real_t *vec = new real_t[N];
    for (len_t i = 0; i < N; i++)
        vec[i] = 0.0;

    // Evaluate equation
    eqn->RebuildEquations(t0, 0, unknowns);
    eqn->Evaluate(uqtyId, vec, unknowns);

    // Store initial value
    this->unknowns->SetInitialValue(uqtyId, vec, t0);

    delete [] vec;
}

/**
 * Initialize the specified unknown quantity by evaluating
 * its initialization function.
 *
 * t0:     Time for which the system should be initialized.
 * uqtyId: ID of the unknown quantity to initialize.
 */
void EqsysInitializer::EvaluateFunction(const real_t t0, const len_t uqtyId) {
    struct initrule *rule = this->rules[uqtyId];
    UnknownQuantityEquation *eqn = this->unknown_equations->at(uqtyId);

    if (!rule->init)
        throw EqsysInitializerException(
            "Unable to initialize '%s': no initialization function given.",
            this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
        );

    // Construct vector for temporary storage of quantity value
    const len_t N = eqn->NumberOfElements();
    real_t *vec = new real_t[N];
    for (len_t i = 0; i < N; i++)
        vec[i] = 0.0;
    
    // Evaluate the unknown quantity
    rule->init(this->unknowns, vec);

    // Store initial value
    this->unknowns->SetInitialValue(uqtyId, vec, t0);

    delete [] vec;
}

