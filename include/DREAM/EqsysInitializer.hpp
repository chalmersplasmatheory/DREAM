#ifndef _DREAM_EQSYS_INITIALIZER_HPP
#define _DREAM_EQSYS_INITIALIZER_HPP

#include <vector>
#include "DREAM/NotImplementedException.hpp"
#include "FVM/config.h"
#include "FVM/FVMException.hpp"
#include "FVM/UnknownQuantityHandler.hpp"

namespace DREAM {
    typedef std::function<void(real_t*)> initfunc_t;

    class EqsysInitializerException : public DREAM::FVM::FVMException {
    public:
        template<typename ... Args>
        EqsysInitializerException(const std::string &msg, Args&& ... args)
            : FVMException(msg, std::forward<Args>(args) ...) {
            AddModule("EqsysInitializer");
        }
    };

    class EqsysInitializer {
    private:
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
        /**
         * An initialization rule which describes how to initialize
         * the specified unknown quantity.
         */
        struct initrule {
            // ID of unknown quantity to which the rule applies
            len_t uqtyId;

            // Rule type
            enum initrule_t type;

            // List of IDs of unknown quantities on which the initial
            // value of this quantity depends
            std::vector<len_t> dependencies;

            // Initialization function (if type = EVAL_FUNCTION)
            initfunc_t *init;

            // DEFAULT PARAMETERS
            // 'true' when this rule has been executed. 'false'
            // until then.
            bool executed = false;
            // 'true' when this rule has been marked as considered
            // for execution. This is used to avoid circular dependencies.
            bool marked = false;
        };

        FVM::UnknownQuantityHandler *unknowns;
        std::unordered_map<len_t, struct initrule> rules;

        std::vector<len_t> ConstructExecutionOrder(struct initrule&);
    public:
        EqsysInitializer(FVM::UnknownQuantityHandler*);
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
        void AddRule(const len_t uqtyId, const enum initrule_t type, initfunc_t *fnc=nullptr);

        template<typename ... Args>
        void AddRule(const len_t uqtyId, const enum initrule_t type, initfunc_t *fnc, Args&& ... deps) {
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

            rules[uqtyId] = {
                uqtyId, type,
                std::vector<len_t>{deps...},
                fnc
            };
        }

        void Execute();
        bool HasRuleFor(const len_t uqtyId) const;
    };
}

#endif/*_DREAM_EQSYS_INITIALIZER_HPP*/
