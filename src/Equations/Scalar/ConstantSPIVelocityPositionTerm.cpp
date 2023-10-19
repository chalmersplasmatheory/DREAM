#include "FVM/config.h"
#include "FVM/Equation/EquationTerm.hpp"
#include "FVM/Equation/PredeterminedParameter.hpp"
#include "FVM/Grid/Grid.hpp"
#include "DREAM/Equations/Scalar/ConstantSPIVelocityPositionTerm.hpp"

using namespace DREAM;

ConstantSPIVelocityPositionTerm::ConstantSPIVelocityPositionTerm(FVM::Grid *g, FVM::UnknownQuantityHandler* u, const real_t *xp,  const real_t *t_delay): FVM::PredeterminedParameter(g){
    id_vp=u->GetUnknownID(OptionConstants::UQTY_V_P);
    len_t id_Yp=u->GetUnknownID(OptionConstants::UQTY_Y_P);
    nShard=u->GetUnknown(id_Yp)->NumberOfMultiples();
    this->currentData = new real_t[3*nShard];
    this->SetData(xp);
    if(t_delay != nullptr)
        this->SetTimeDelay(t_delay);
}
ConstantSPIVelocityPositionTerm::~ConstantSPIVelocityPositionTerm(){
    delete [] currentData;
    delete [] initial_data;
    delete [] t_delay;
}


void ConstantSPIVelocityPositionTerm::SetData(const real_t* xp_init, bool copy){
    if(copy){
        real_t *initial_data_xp=new real_t[3*this->nShard];
        for(len_t ip=0;ip<nShard;ip++){
            initial_data_xp[3*ip]=xp_init[3*ip];
            initial_data_xp[3*ip+1]=xp_init[3*ip+1];
            initial_data_xp[3*ip+2]=xp_init[3*ip+2];
        }
        this->initial_data=initial_data_xp;
    }else{
        this->initial_data=xp_init;
    }
}
void ConstantSPIVelocityPositionTerm::SetTimeDelay(const real_t *t_delay, bool copy){
    if(copy){
        real_t *data_tmp=new real_t[nShard];
        for(len_t ip=0;ip<nShard;ip++){
            data_tmp[ip]=t_delay[ip];
        }
        this->t_delay=data_tmp;
    }else{
        this->t_delay=t_delay;
    }
}
void ConstantSPIVelocityPositionTerm::Rebuild(const real_t t, const real_t dt, FVM::UnknownQuantityHandler* u){
    vp=u->GetUnknownData(id_vp);
    for(len_t ip=0;ip<nShard;ip++){
        if(t+dt>t_delay[ip]){
            currentData[3*ip]=initial_data[3*ip]+vp[3*ip]*(t+dt-t_delay[ip]);
            currentData[3*ip+1]=initial_data[3*ip+1]+vp[3*ip+1]*(t+dt-t_delay[ip]);
            currentData[3*ip+2]=initial_data[3*ip+2]+vp[3*ip+2]*(t+dt-t_delay[ip]);
        }else{
            currentData[3*ip]=initial_data[3*ip];
            currentData[3*ip+1]=initial_data[3*ip+1];
            currentData[3*ip+2]=initial_data[3*ip+2];
        }
    }
}
