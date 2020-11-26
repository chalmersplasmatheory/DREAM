#include "DREAM/Equations/Fluid/MaxwellianCollisionalEnergyTransferTerm.hpp"
#include "DREAM/Constants.hpp"

/**
 * Implementation of an equation terms which represents the rate of energy transfer (in SI)
 * between two maxwellians of densities n_i, n_j and heat content W_i, W_j with W = 1.5nT
 * representing particles with masses m_i, m_j and charge numbers Z_i and Z_j.
 * 
 * It can describe both ion-ion and ion-electron heat transfer, only differing in 
 * definition by the Coulomb logarithm used, controlled with the `bool isEI` flag.
 * 
 * One instance of this class represents the interaction between just two species.
 */

using namespace DREAM;

/**
 * Constructor.
 * The unknown id's defines which quantities are targets and which are "incident" 
 * (should be either ions or electrons). The offsets should be 0 for electrons and
 * the ion index for ions. 
 */
MaxwellianCollisionalEnergyTransferTerm::MaxwellianCollisionalEnergyTransferTerm(
            FVM::Grid *g, 
            const len_t id_ni, const len_t id_Wi, const len_t Zi, const real_t mi, const len_t offset_i, 
            const len_t id_nj, const len_t id_Wj, const len_t Zj, const real_t mj, const len_t offset_j, 
            FVM::UnknownQuantityHandler *u, CoulombLogarithm *lnL, bool isEI, real_t scaleFactor
) : FVM::EquationTerm(g), id_ni(id_ni), id_nj(id_nj), id_Wi(id_Wi), id_Wj(id_Wj),
    Zi(Zi), Zj(Zj), offset_i(offset_i), offset_j(offset_j), mi(mi), mj(mj), 
    unknowns(u), lnLambda(lnL), isEI(isEI)
{
    real_t ec = Constants::ec;
    real_t eps0 = Constants::eps0;
    real_t ec2eps = ec*ec / eps0;
    this->constPreFactor = scaleFactor*sqrt(3.0/M_PI) * Zi*Zi*Zj*Zj 
                            *ec2eps*ec2eps*sqrt(mi*mj) / (4*M_PI);

    this->lnLambda_settings = new CollisionQuantity::collqty_settings;
    this->lnLambda_settings->lnL_type = isEI ? 
        OptionConstants::COLLQTY_LNLAMBDA_THERMAL : OptionConstants::COLLQTY_LNLAMBDA_ION_ION;

    this->id_ions  = u->GetUnknownID(OptionConstants::UQTY_ION_SPECIES);
    this->id_Tcold = u->GetUnknownID(OptionConstants::UQTY_T_COLD);
}

/**
 * Sets the jacobian block of this equation term
 */
void MaxwellianCollisionalEnergyTransferTerm::SetJacobianBlock(const len_t /*uqtyId*/, const len_t derivId, FVM::Matrix *jac, const real_t*){
    if(derivId != id_ni && derivId != id_nj && derivId != id_Wi && derivId != id_Wj)
        return;

    const real_t *ni = unknowns->GetUnknownData(id_ni)+offset_i*nr;
    const real_t *nj = unknowns->GetUnknownData(id_nj)+offset_j*nr;
    const real_t *Wi = unknowns->GetUnknownData(id_Wi)+offset_i*nr;
    const real_t *Wj = unknowns->GetUnknownData(id_Wj)+offset_j*nr;
    const real_t *lnL = isEI ? lnLambda->GetLnLambdaT() : lnLambda->GetLnLambdaII();
    
    for(len_t ir=0; ir<nr; ir++){
        real_t njWi = nj[ir]*Wi[ir];
        real_t niWj = ni[ir]*Wj[ir];
        if(ni[ir]==0 || nj[ir]==0 || niWj==njWi) // last one (Ti=Tj) deals with Wi=Wj=0 which is singular below
            continue;
        real_t pre = sqrt(ni[ir]*nj[ir])*ni[ir]*nj[ir];
        real_t up = niWj - njWi;
        real_t down = mj*njWi + mi*niWj;
        real_t vec = lnL[ir] * constPreFactor * pre * up / (down*sqrt(down));

        len_t ii = offset_i*nr + ir;
        len_t jj = offset_j*nr + ir;
        if(derivId==id_ni)
            jac->SetElement(ii,ii, vec*( 1.5/ni[ir] + Wj[ir]*( 1.0/up - 1.5*mi/down) ) );
        if(derivId==id_nj) 
            jac->SetElement(ii,jj, vec*( 1.5/nj[ir] + Wi[ir]*(-1.0/up - 1.5*mj/down) ) );
        if(derivId==id_Wi)
            jac->SetElement(ii,ii, vec*nj[ir]*(-1.0/up - 1.5*mj/down ) );
        if(derivId==id_Wj) 
            jac->SetElement(ii,jj, vec*ni[ir]*( 1.0/up - 1.5*mi/down ) );
        // Below: lnLambda derivatives
        if(derivId==id_Tcold || derivId==id_ions)
            for(len_t n=0; n<unknowns->GetUnknown(derivId)->NumberOfMultiples(); n++)
                jac->SetElement(ii,n*nr+ir, vec/lnL[ir]*lnLambda->evaluatePartialAtP(ir,0,derivId,n,lnLambda_settings) );
    }
}

/**
 * Set matrix elements of this equation term; this equation term is treated explicitly
 * and everything is put in the right-hand side, as there is no clear unknown that
 * represents this equation term especially well.
 */
void MaxwellianCollisionalEnergyTransferTerm::SetMatrixElements(FVM::Matrix* /*mat*/, real_t *rhs){
    SetVectorElements(rhs, nullptr);
}

/**
 * Set vector elements of this equation term. 
 */
void MaxwellianCollisionalEnergyTransferTerm::SetVectorElements(real_t *vec, const real_t*){
    const real_t *ni = unknowns->GetUnknownData(id_ni)+offset_i*nr;
    const real_t *nj = unknowns->GetUnknownData(id_nj)+offset_j*nr;
    const real_t *Wi = unknowns->GetUnknownData(id_Wi)+offset_i*nr;
    const real_t *Wj = unknowns->GetUnknownData(id_Wj)+offset_j*nr;

    const real_t *lnL = isEI ? lnLambda->GetLnLambdaT() : lnLambda->GetLnLambdaII();
    for(len_t ir=0; ir<nr; ir++){
        real_t njWi = nj[ir]*Wi[ir];
        real_t niWj = ni[ir]*Wj[ir];
        if(ni[ir]==0 || nj[ir]==0 || niWj==njWi)
            continue;
        real_t pre = sqrt(ni[ir]*nj[ir])*ni[ir]*nj[ir];
        real_t up = niWj - njWi;
        real_t down = mj*njWi + mi*niWj;
        vec[offset_i*nr + ir] += lnL[ir] * constPreFactor * pre * up / (down*sqrt(down));
    }
}
