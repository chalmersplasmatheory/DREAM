#include "DREAM/Equations/Fluid/IonisationHeatingTerm.hpp"

/**
 * Implementation of 
 */


using namespace DREAM;


IonisationHeatingTerm::IonisationHeatingTerm(
    FVM::Grid* g, FVM::UnknownQuantityHandler *u, 
    IonHandler *ionHandler, ADAS *adas, NIST *nist) 
    : FVM::DiagonalComplexTerm(g,u) 
{
    SetName("IonisationHeatingTerm");

    this->adas = adas;
    this->ionHandler = ionHandler;
    this->nist = nist;

    this->id_ncold = unknowns->GetUnknownID(OptionConstants::UQTY_N_COLD);
    this->id_nhot  = unknowns->GetUnknownID(OptionConstants::UQTY_N_HOT);
    this->id_Tcold = unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD);
    this->id_ni    = unknowns->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);

    AddUnknownForJacobian(unknowns, id_ncold);
    AddUnknownForJacobian(unknowns, id_nhot);
    AddUnknownForJacobian(unknowns, id_ni);
    
}


void IonisationHeatingTerm::AllocateInitParams()
{
    nr_init = nr;
    if(T_cold_init != nullptr)
        delete [] T_cold_init;
    T_cold_init = new real_t[nr];
/*
    if(n_cold_init != nullptr)
        delete [] n_cold_init;
    n_cold_init = new real_t[nr];
  */  
}


void IonisationHeatingTerm::SetWeights() 
{
    len_t NCells = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();
    
    real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *n_hot  = unknowns->GetUnknownData(id_nhot);
    real_t *T_cold = unknowns->GetUnknownData(id_Tcold);
    real_t *n_i    = unknowns->GetUnknownData(id_ni);
    
    if(nr_init != nr){
        AllocateInitParams();
        for(len_t ir=0; ir<nr; ir++)
            T_cold_init[ir] = T_cold[ir];
    }

    for (len_t i = 0; i < NCells; i++)
            weights[i] = 0;


    for(len_t iz = 0; iz<nZ; iz++){
        ADASRateInterpolator *scd = adas->GetSCD(Zs[iz]);        
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
            len_t nMultiple = ionHandler->GetIndex(iz,Z0);
            for (len_t i = 0; i < NCells; i++){
                real_t Ion =  scd->Eval(Z0, n_cold[i]+n_hot[i], T_cold_init[i]);
                real_t I_energy = nist->GetIonizationEnergy(Zs[iz],Z0) * Constants::ec;
                real_t ni = n_i[nMultiple*NCells + i];
                weights[i] += -ni * I_energy * Ion;
            }
        }
    }
}

void IonisationHeatingTerm::SetDiffWeights(len_t derivId, len_t /*nMultiples*/){
    len_t NCells = grid->GetNCells();
    len_t nZ = ionHandler->GetNZ();
    const len_t *Zs = ionHandler->GetZs();

    real_t *n_cold = unknowns->GetUnknownData(id_ncold);
    real_t *n_hot  = unknowns->GetUnknownData(id_nhot);
    real_t *n_i    = unknowns->GetUnknownData(id_ni);

    if(derivId == id_ni)
        for(len_t iz = 0; iz<nZ; iz++){
            ADASRateInterpolator *scd = adas->GetSCD(Zs[iz]);        
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                len_t nMultiple = ionHandler->GetIndex(iz,Z0);
                for (len_t i = 0; i < NCells; i++){
                    real_t Ion =  scd->Eval(Z0, n_cold[i]+n_hot[i], T_cold_init[i]);
                    real_t I_energy = nist->GetIonizationEnergy(Zs[iz],Z0) * Constants::ec;
                    diffWeights[NCells*nMultiple + i] = -I_energy * Ion;
                }
            }
        }
    else if( (derivId == id_ncold) || (derivId == id_nhot))
        for(len_t iz = 0; iz<nZ; iz++){
            ADASRateInterpolator *scd = adas->GetSCD(Zs[iz]);        
            for(len_t Z0 = 0; Z0<=Zs[iz]; Z0++){
                len_t nMultiple = ionHandler->GetIndex(iz,Z0);
                for (len_t i = 0; i < NCells; i++){
                    real_t dIon =  scd->Eval_deriv_n(Z0, n_cold[i]+n_hot[i], T_cold_init[i]);
                    real_t I_energy = nist->GetIonizationEnergy(Zs[iz],Z0) * Constants::ec;
                    real_t ni = n_i[nMultiple*NCells + i];
                    diffWeights[i] = -ni * I_energy * dIon;
                }
            }
        }
}

