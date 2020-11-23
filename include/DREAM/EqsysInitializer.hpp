#ifndef _DREAM_EQSYS_INITIALIZER_HPP
#define _DREAM_EQSYS_INITIALIZER_HPP

#include <string>
#include <vector>
#include <softlib/SFile.h>
#include "DREAM/Equations/CollisionQuantityHandler.hpp"
#include "DREAM/IonHandler.hpp"
#include "DREAM/NotImplementedException.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
#include "DREAM/UnknownQuantityEquation.hpp"
#include "FVM/config.h"
#include "FVM/FVMException.hpp"
#include "FVM/Interpolator3D.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    typedef std::function<void(FVM::UnknownQuantityHandler*, real_t*)> initfunc_t;

    class EqsysInitializerException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        EqsysInitializerException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("EqsysInitializer");
        }
    };

    class EqsysInitializer {
    public:
        enum initrule_t {
            // Evaluated from equation for unknown (which is of the form
            //   c*x - f(y) = 0
            // where 'x' denotes the unknown quantity, c is (possibly radially
            // varying) constant, and 'y' denotes any number of _other_ unknown
            // quantities (i.e. NOT x).
            INITRULE_EVAL_EQUATION,

            // Evaluated using a function which is specified when creating
            // the rule.
            INITRULE_EVAL_FUNCTION,

            // Evaluated by solving a system of equations, setting any
            // transient terms to zero.
            INITRULE_STEADY_STATE_SOLVE
        };

        // Special types
        static const int_t
            COLLQTYHDL_HOTTAIL,
            COLLQTYHDL_RUNAWAY,
            RUNAWAY_FLUID;

    private:
        /**
         * An initialization rule which describes how to initialize
         * the specified unknown quantity.
         */
        struct initrule {
            initrule() {}
            initrule(int_t u, enum initrule_t t, const std::vector<int_t>& l, initfunc_t f)
                : uqtyId(u), type(t), dependencies(l), init(f) { }
            initrule(int_t u, enum initrule_t t, const std::initializer_list<int_t>& l, initfunc_t f)
                : uqtyId(u), type(t), dependencies(l), init(f) { }

            // ID of unknown quantity to which the rule applies
            int_t uqtyId;

            // Rule type
            enum initrule_t type;

            // List of IDs of unknown quantities on which the initial
            // value of this quantity depends
            std::vector<int_t> dependencies;

            // Initialization function (if type = EVAL_FUNCTION)
            initfunc_t init;

            // DEFAULT PARAMETERS
            // 'true' when this rule has been executed. 'false'
            // until then.
            bool executed = false;
            // 'true' when this rule has been marked as considered
            // for execution. This is used to avoid circular dependencies.
            bool marked = false;
        };

        FVM::UnknownQuantityHandler *unknowns;
        std::vector<UnknownQuantityEquation*> *unknown_equations;
        std::unordered_map<int_t, struct initrule*> rules;

        std::vector<int_t> ConstructExecutionOrder(struct initrule*);

        FVM::Grid *fluidGrid, *hottailGrid, *runawayGrid;
        enum OptionConstants::momentumgrid_type
            hottail_type, runaway_type;

        CollisionQuantityHandler *cqhHottail=nullptr, *cqhRunaway=nullptr;
        IonHandler *ionHandler=nullptr;
        RunawayFluid *runawayFluid=nullptr;

        void __InitTR(
            FVM::UnknownQuantity*, const real_t, const int_t,
            const len_t, const real_t*, const real_t*, const sfilesize_t*
        );
        void __InitTRmult(
            FVM::UnknownQuantity*, const real_t, const int_t,
            const len_t, const real_t*, const real_t*,
            const sfilesize_t*
        );
        void __InitTR2P(
            FVM::UnknownQuantity*, const real_t, const int_t,
            const real_t*, const real_t*, const real_t*, const real_t*,
            const sfilesize_t*, enum OptionConstants::momentumgrid_type,
            enum OptionConstants::momentumgrid_type
        );
    public:
        EqsysInitializer(
            FVM::UnknownQuantityHandler*, std::vector<UnknownQuantityEquation*>*,
            FVM::Grid*, FVM::Grid*, FVM::Grid*,
            enum OptionConstants::momentumgrid_type,
            enum OptionConstants::momentumgrid_type
        );
        ~EqsysInitializer();

        /**
         * Add a new initialization rule.
         *
         * uqtyId: ID of the unknown quantity to which the rule applies.
         * type:   Rule type (see 'enum initrule_t').
         * fnc:    Optional initialization function (which is called when
         *         'type = INITRULE_EVAL_FUNCTION' instead of evaluating the
         *         unknown's equation).
         * deps:   Optional list of IDs of unknown quantities on which this
         *         initialization rule depends (i.e. a list of quantities
         *         which must be initialized before initializing 'uqtyId').
         */
        void AddRule(const std::string& uqtyName, const enum initrule_t type, initfunc_t fnc=nullptr);
        void AddRule(const int_t uqtyId, const enum initrule_t type, initfunc_t fnc=nullptr);

        template<typename ... Args>
        void AddRule(const std::string& uqtyName, const enum initrule_t type, initfunc_t fnc, Args&& ... deps) {
            const int_t uqtyId = (int_t)this->unknowns->GetUnknownID(uqtyName);
            this->AddRule(uqtyId, type, fnc, deps...);
        }
        template<typename ... Args>
        void AddRule(const int_t uqtyId, const enum initrule_t type, initfunc_t fnc, Args&& ... deps) {
            if (this->HasRuleFor(uqtyId))
                throw EqsysInitializerException(
                    "An initialization rule for '%s' has already been created.",
                    this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
                );
            if (type == INITRULE_STEADY_STATE_SOLVE)
                throw NotImplementedException(
                    "%s: No support for solving for initial steady state implemented yet.",
                    this->unknowns->GetUnknown(uqtyId)->GetName().c_str()
                );

            rules[uqtyId] = new struct initrule(
                uqtyId, type,
                std::vector<int_t>{static_cast<int_t>(deps)...},
                fnc
            );
        }
        void RemoveRule(const int_t);

        void EvaluateEquation(const real_t, const int_t);
        void EvaluateFunction(const real_t, const int_t);

        void Execute(const real_t);
        bool HasRuleFor(const int_t uqtyId) const;
        void InitializeFromOutput(const std::string&, const real_t, int_t, IonHandler*, std::vector<std::string>&);
        void InitializeFromOutput(SFile*, const real_t, int_t, IonHandler*, std::vector<std::string>&);
        void VerifyAllInitialized() const;

        void SetHottailCollisionHandler(CollisionQuantityHandler *cqh)
        { this->cqhHottail = cqh; }
        void SetRunawayCollisionHandler(CollisionQuantityHandler *cqh)
        { this->cqhRunaway = cqh; }
        void SetRunawayFluid(RunawayFluid *ref)
        { this->runawayFluid = ref; }
        void SetIonHandler(IonHandler *ih)
        { this->ionHandler = ih; }
    };
}

#endif/*_DREAM_EQSYS_INITIALIZER_HPP*/
