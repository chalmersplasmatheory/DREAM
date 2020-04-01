#ifndef _DREAM_EQUATIONS_COLLISION_FREQUENCY_CREATOR_HPP
#define _DREAM_EQUATIONS_COLLISION_FREQUENCY_CREATOR_HPP

#include "FVM/config.h"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Constants.hpp"
#include <cmath>
namespace DREAM {
    class CollisionFrequencyCreator{
    
    protected:
        

        real_t *coulombLog_c = nullptr;
        
        FVM::Grid *grid;
        EquationSystem *eqSys;
        
    public:
        CollisionFrequencyCreator(FVM::Grid*,EquationSystem*);
        ~CollisionFrequencyCreator();

        const real_t constPreFactor = 4*Constants::pi 
                           *Constants::r0*Constants::r0
                           *Constants::c;
        
        void RebuildLnLambda_c(real_t*, real_t*);

        real_t *n_cold;
        real_t *T_cold;
        real_t Zeff=1; //todo

        real_t nu_s(len_t, real_t);
        real_t nu_D(len_t, real_t);
        real_t lnLambda_ee(len_t, real_t);
        real_t lnLambda_ei(len_t, real_t);
        
        real_t lnLambda_c(len_t ir){
            return coulombLog_c[ir]; }

        virtual void Refresh();

    };
}

#endif/*_DREAM_EQUATIONS_COLLISION_FREQUENCY_CREATOR_HPP*/
