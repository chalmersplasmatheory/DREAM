/**
 * Implementation of the "CollisionFrequencyCreator" class, which 
 * manages the calculation of collision frequencies such as
 * nu_s, nu_D and lnLambda, based on data from the unknown_variables of an
 * EquationSystem. 
 * In the future, it could also load atomic physics data and 
 * generate ionization rates and such  
 */

//#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/EquationSystem.hpp"
#include "DREAM/Equations/CollisionFrequencyCreator.hpp"
#include "cmath"
//#include <iostream>
using namespace DREAM;

/**
 * Constructor.
 */
CollisionFrequencyCreator::CollisionFrequencyCreator(FVM::Grid *g,EquationSystem *es) {
    this->coulombLog_c = new real_t[g->GetNr()];
    this->grid = g;
    this->eqSys = es;
}


/**
 * Destructor.
 */
CollisionFrequencyCreator::~CollisionFrequencyCreator(){
    delete [] coulombLog_c;
}

void CollisionFrequencyCreator::RebuildLnLambda_c(real_t *n, real_t *T){
    len_t nr = this->grid->GetNr();
    for (len_t ir = 0; ir<nr; ir++)
        this->coulombLog_c[ir] = 14.6 + 0.5*log(T[ir]/(n[ir]/1e20));
}

real_t CollisionFrequencyCreator::lnLambda_ee(len_t ir, real_t p){
    real_t gamma = sqrt(1+p*p);
    return CollisionFrequencyCreator::lnLambda_c(ir) + log(sqrt(gamma-1)); 
}
real_t CollisionFrequencyCreator::lnLambda_ei(len_t ir, real_t p){
    return CollisionFrequencyCreator::lnLambda_c(ir) + log(sqrt(2)*p); 
}
        
        

//Calculates slowing-down frequency. Later this will load 
//atomic data from Linnea's tables and sum over ion species
real_t CollisionFrequencyCreator::nu_s(len_t ir, real_t p){
    return this->constPreFactor*this->n_cold[ir]*CollisionFrequencyCreator::lnLambda_ee(ir, p) 
                    * (1+p*p)/(p*p*p); 
}


//Calculates pitch-angle scattering frequency. Later this will  
//load atomic data from Linnea's tables and sum over ion species
real_t CollisionFrequencyCreator::nu_D(len_t ir, real_t p){
    real_t gamma = sqrt(1+p*p);
    return this->constPreFactor * (1+this->Zeff) * this->n_cold[ir]
                    *CollisionFrequencyCreator::lnLambda_ei(ir, p) * gamma/(p*p*p);
}

void CollisionFrequencyCreator::Refresh(){
    len_t id_ncold = eqSys->GetUnknownID("n_cold");
    len_t id_Tcold = eqSys->GetUnknownID("T_cold");
    this->n_cold   = eqSys->GetUnknownData(id_ncold);
    this->T_cold   = eqSys->GetUnknownData(id_Tcold);
    //treat impurity densities, charge states to find Zeff etc
    //this->Zeff = ...
    CollisionFrequencyCreator::RebuildLnLambda_c(this->n_cold,this->T_cold); 
}

        
