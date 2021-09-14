#include "DREAM/Equations/Fluid/CollisionalEnergyTransferKineticTerm.hpp"

using namespace DREAM;


/**
 * Constructor. 
 */
CollisionalEnergyTransferKineticTerm::CollisionalEnergyTransferKineticTerm(
    FVM::Grid *densityGrid, FVM::Grid *distributionGrid, len_t id_n, len_t id_f,
    CollisionQuantityHandler* cqh, FVM::UnknownQuantityHandler *u, 
    enum OptionConstants::momentumgrid_type mgt, real_t sf,
    real_t pThreshold, pThresholdMode pMode
) : MomentQuantity(densityGrid, distributionGrid, id_n, id_f,u, 
    pThreshold, pMode), collQtyHandler(cqh), mgtype(mgt), scaleFactor(sf) {

    SetName("CollisionalEnergyTransferKineticTerm");

    /**
     * Using "SUPERTHERMAL" collision setting to evaluate energy transfer, 
     * and completely screened type (so that only colliding with n_cold).
     * Neglecting bremsstrahlung and using energy-dependent lnLambda.
     */
    this->collQtySetting = new CollisionQuantity::collqty_settings;
    this->collQtySetting->collfreq_mode = OptionConstants::COLLQTY_COLLISION_FREQUENCY_MODE_SUPERTHERMAL;
    this->collQtySetting->collfreq_type = OptionConstants::COLLQTY_COLLISION_FREQUENCY_TYPE_COMPLETELY_SCREENED;
    this->collQtySetting->lnL_type      = OptionConstants::COLLQTY_LNLAMBDA_ENERGY_DEPENDENT;    

    AddUnknownForJacobian(unknowns, OptionConstants::UQTY_T_COLD);
    AddUnknownForJacobian(unknowns, OptionConstants::UQTY_N_COLD);
    AddUnknownForJacobian(unknowns, OptionConstants::UQTY_ION_SPECIES);
    AllocateDiffIntegrand();
}

CollisionalEnergyTransferKineticTerm::~CollisionalEnergyTransferKineticTerm(){
    delete [] collQtySetting;
}


/**
 * The collisional energy transfer density is approximately set by p*v*nu_s*f,
 * where nu_s is the collisional friction against free electrons only
 */
void CollisionalEnergyTransferKineticTerm::Rebuild(const real_t, const real_t, FVM::UnknownQuantityHandler*){
    real_t p,v,nu_s;
    len_t offset = 0;
    for(len_t ir = 0; ir<fGrid->GetNr(); ir++){
        FVM::MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
        real_t np1 = mg->GetNp1();
        real_t np2 = mg->GetNp2();
        if(mgtype == OptionConstants::MOMENTUMGRID_TYPE_PXI){
            for(len_t i = 0; i<np1; i++){
                p = mg->GetP1(i);
                v = Constants::c * p/sqrt(1+p*p);
                nu_s = collQtyHandler->GetNuS()->evaluateAtP(ir,p,collQtySetting);
                real_t toSet = scaleFactor * Constants::me * Constants::c *v*p*nu_s;
                for(len_t j = 0; j<np2; j++){
                    len_t ind = offset + j*np1 + i;
                    this->integrand[ind] = toSet;
                }
            }
        } else {
            for(len_t ip1 = 0; ip1<np1; ip1++)
                for(len_t ip2 = 0; ip2<np2; ip2++){
                    len_t ind = offset + ip2*np1 + ip1;
                    p = mg->GetP(ip1,ip2);
                    v = Constants::c * p/mg->GetGamma(ip1,ip2);
                    nu_s = collQtyHandler->GetNuS()->evaluateAtP(ir,p,collQtySetting);
                    this->integrand[ind] = scaleFactor * Constants::me * Constants::c *v*p*nu_s;
                }
        }
        offset += np1*np2;
    }
}

void CollisionalEnergyTransferKineticTerm::SetDiffIntegrand(len_t derivId){
    real_t p,v,dNu_s;
    len_t nMultiples = unknowns->GetUnknown(derivId)->NumberOfMultiples();

    len_t offset = 0;
    for(len_t n=0; n<nMultiples; n++)
        for(len_t ir = 0; ir<fGrid->GetNr(); ir++){
            FVM::MomentumGrid *mg = fGrid->GetMomentumGrid(ir);
            real_t np1 = mg->GetNp1();
            real_t np2 = mg->GetNp2();

            if(mgtype == OptionConstants::MOMENTUMGRID_TYPE_PXI){
                for(len_t i = 0; i<np1; i++){
                    p = mg->GetP1(i);
                    v = Constants::c * p/sqrt(1+p*p);
                    dNu_s = collQtyHandler->GetNuS()->evaluatePartialAtP(ir,p,derivId,n,collQtySetting);
                    real_t toAdd = scaleFactor * Constants::me * Constants::c *v*p*dNu_s;
                    for(len_t j = 0; j<np2; j++){ // Possible optimization: for pxi-grid, evaluate dNuS outside of i2 loop
                        len_t ind = offset + j*np1 + i;
                        this->diffIntegrand[ind] += toAdd;
                    }
                }
            } else{
                for(len_t ip1 = 0; ip1<np1; ip1++)
                    for(len_t ip2 = 0; ip2<np2; ip2++){ // Possible optimization: for pxi-grid, evaluate dNuS outside of i2 loop
                        len_t ind = offset + ip2*np1 + ip1;
                        p = mg->GetP(ip1,ip2);
                        v = Constants::c * p/mg->GetGamma(ip1,ip2);
                        dNu_s = collQtyHandler->GetNuS()->evaluatePartialAtP(ir,p,derivId,n,collQtySetting);
                        this->diffIntegrand[ind] += scaleFactor * Constants::me * Constants::c *v*p*dNu_s;
                    }
            }
            offset += np1*np2;
        }
}
