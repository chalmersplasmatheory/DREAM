#ifndef _DREAM_EQUATIONS_KINETIC_BOLTZMANN_KNOCK_ON_OPERATOR_HPP
#define _DREAM_EQUATIONS_KINETIC_BOLTZMANN_KNOCK_ON_OPERATOR_HPP

#include "FVM/Equation/EquationTerm.hpp"
#include "DREAM/Settings/OptionConstants.hpp"
namespace DREAM {
    class BoltzmannKnockOnOperator
        : public FVM::EquationTerm {

    private:
        enum pitch_interval {
            PITCH_INTERVAL_LOWER,
            PITCH_INTERVAL_MIDDLE,
            PITCH_INTERVAL_UPPER
        };
        struct ParametersForPitchTerm {
            real_t xi0Prime;
            real_t xiStar;
            real_t thetaStar;
            real_t xi0_jp;
            real_t xi0_jm;
            len_t ir;
            FVM::RadialGrid *rGrid;
        };
        FVM::UnknownQuantityHandler *unknowns;
        len_t id_ntot;
        real_t ***SourceMatrix = nullptr;

        real_t EvaluatePitchIntegrandAtTheta(real_t theta, void *par);
        real_t PitchIntegrandAtXi(
            real_t xi0_lo, real_t xi0_up, real_t xiPrime, real_t thetaPrime, real_t xiStar, real_t thetaStar
        );
        real_t PitchIntegrandAtanFunc(
            real_t xi_lower, real_t xi_upper, real_t xiPrime, real_t xiStar,
            real_t theta_lower, real_t theta_upper, real_t thetaPrime, real_t thetaStar,
            pitch_interval interval
        );
        void Deallocate();
        void Allocate();
    protected:

    public:
        BoltzmannKnockOnOperator(FVM::Grid*, FVM::UnknownQuantityHandler*);
        ~BoltzmannKnockOnOperator();

        virtual void SetMatrixElements(FVM::Matrix*, real_t*) override;
        virtual void SetVectorElements(real_t*, const real_t*) override;
        virtual void SetJacobianBlock(const len_t uqtyId, const len_t derivId, FVM::Matrix *jac, const real_t* x) override;

        virtual len_t GetNumberOfNonZerosPerRow() const override;

        virtual void Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*) override;
        virtual bool GridRebuilt() override;

        real_t EvaluateMollerDifferentialCrossSection(real_t gamma, real_t gamma1);
        real_t EvaluateIntegratedMollerDCS(real_t gamma_l, real_t gamma_u, real_t gamma1);
    };
}

#endif/*_DREAM_EQUATIONS_KINETIC_BOLTZMANN_KNOCK_ON_OPERATOR_HPP*/


