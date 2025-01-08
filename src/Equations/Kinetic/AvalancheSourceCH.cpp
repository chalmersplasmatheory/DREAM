/**
 * Implementation of the Chiu-Harvey avalanche
 * source term, which takes the quadratic form
 *     T = S(r,p) * n_tot(r,t) * f_re(r,t),
 * where S is only a function of phase-space coordinates.
 */

#include "DREAM/Equations/Kinetic/AvalancheSourceCH.hpp"
#include "DREAM/Constants.hpp"

using namespace DREAM;

/**
 * Constructor.
 */
AvalancheSourceCH::AvalancheSourceCH(
    FVM::Grid *kineticGrid, FVM::UnknownQuantityHandler *u,
    real_t pCutoff, real_t scaleFactor, CHSourceMode sm, 
    bool isRunawayGrid, FVM::Grid* runawayGrid
) : FluidSourceTerm(kineticGrid, u), scaleFactor(scaleFactor),
    sourceMode(sm), isRunawayGrid(isRunawayGrid), runawayGrid(runawayGrid)
{
    SetName("AvalancheSourceCH");
    
    id_Efield = unknowns->GetUnknownID(OptionConstants::UQTY_E_FIELD);
    this->pCutoff = pCutoff;
    
    real_t e = Constants::ec;
    real_t epsmc = Constants::eps0 * Constants::me * Constants::c;
    this->preFactor = (e*e*e*e)/(2*M_PI*epsmc*epsmc*Constants::c);

    if (isRunawayGrid) {
        id_fgrid = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);
        AddUnknownForJacobian(u, id_fgrid);
    } else if (runawayGrid != nullptr) {
        htgridWithREgrid = true;
        id_fre = unknowns->GetUnknownID(OptionConstants::UQTY_F_RE);  
        AddUnknownForJacobian(u, id_fre);  
        id_fgrid = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
        AddUnknownForJacobian(u, id_fgrid);
    } else {
        id_fgrid = unknowns->GetUnknownID(OptionConstants::UQTY_F_HOT);
        AddUnknownForJacobian(u, id_fgrid);
    }

    this->pIn_Indices = new int_t[this->grid->GetNCells()];
    this->f_pIndices = new std::vector<len_t>[nr*this->grid->GetMomentumGrid(0)->GetNp1()];
    this->f_xiIndices = new std::vector<len_t>[nr*this->grid->GetMomentumGrid(0)->GetNp1()];
    if (htgridWithREgrid){
        this->f_re_pIndices = new std::vector<len_t>[nr*this->runawayGrid->GetMomentumGrid(0)->GetNp1()];
        this->f_re_xiIndices = new std::vector<len_t>[nr*this->runawayGrid->GetMomentumGrid(0)->GetNp1()];
    }
    SetPinIndices();
}

AvalancheSourceCH::~AvalancheSourceCH(){
    delete [] this->pIn_Indices;
    delete [] this->f_pIndices;
    delete [] this->f_xiIndices;
    if (htgridWithREgrid){
        delete [] this->f_re_pIndices;
        delete [] this->f_re_xiIndices;
    }
}

void AvalancheSourceCH::SetPinIndices(){
    len_t offset = 0;
    for(len_t ir=0; ir<nr; ir++){
        for(len_t i=0; i<n1[ir]; i++)
            for(len_t j=0; j<n2[ir]; j++){
                len_t ind = offset + n1[ir]*j + i;
                
                real_t p = this->grid->GetMomentumGrid(ir)->GetP1(i);
                real_t xi = abs(this->grid->GetMomentumGrid(ir)->GetP2(j));
                real_t gamma = sqrt(p*p + 1.);
                real_t p_in = 2 * xi * p / (xi*xi * (gamma + 1) - gamma + 1);
                
                int_t i_in = -1;
                if (htgridWithREgrid 
                        && p_in > this->grid->GetMomentumGrid(ir)->GetP1_f(n1[ir]) 
                        && p_in < this->runawayGrid->GetMomentumGrid(ir)->GetP1_f(this->runawayGrid->GetMomentumGrid(ir)->GetNp1())){ 
                    
                    i_in = 0;
                    for(len_t i_temp=1; i_temp<this->runawayGrid->GetMomentumGrid(ir)->GetNp1(); i_temp++){
                        if (p_in > this->runawayGrid->GetMomentumGrid(ir)->GetP1_f(i_temp) ){
                            i_in = i_temp;
                        }
                    }
                    
                    this->pIn_Indices[ind] = i_in + n1[ir];
                    this->f_re_pIndices[i_in].push_back(i);
                    this->f_re_xiIndices[i_in].push_back(j);
                } else if (p_in < this->grid->GetMomentumGrid(ir)->GetP1_f(n1[ir])) {

                    i_in = 0;
                    for(len_t i_temp=1; i_temp<n1[ir]; i_temp++){
                        if (p_in > this->grid->GetMomentumGrid(ir)->GetP1_f(i_temp) ){
                            i_in = i_temp;
                        }
                    }
                    
                    this->pIn_Indices[ind] = i_in;
                    this->f_pIndices[i_in].push_back(i);
                    this->f_xiIndices[i_in].push_back(j);
                } else {
                    this->pIn_Indices[ind] = i_in;
                }
            }
        offset += n1[ir]*n2[ir];
    }

}

/**
 * Evaluates the constant (only grid dependent) source-shape function S(r,p)
 */
real_t AvalancheSourceCH::EvaluateCHSource(len_t ir, len_t i, len_t j){
    real_t pm = this->grid->GetMomentumGrid(ir)->GetP1_f(i);
    real_t pp = this->grid->GetMomentumGrid(ir)->GetP1_f(i+1);
    
    real_t cutFactor = 1.;
    if(pp<=pCutoff)
        return 0;
    else if(pm<pCutoff)
        cutFactor = (pp - pCutoff)/(pp - pm); 
    
    real_t p     = this->grid->GetMomentumGrid(ir)->GetP1(i);
    real_t xi    = this->grid->GetMomentumGrid(ir)->GetP2(j);
    real_t gamma = sqrt(p*p + 1);
    real_t pMax  = this->grid->GetMomentumGrid(ir)->GetP1(grid->GetMomentumGrid(ir)->GetNp1()-1);
    if (htgridWithREgrid){
        pMax = this->runawayGrid->GetMomentumGrid(ir)->GetP1(runawayGrid->GetMomentumGrid(ir)->GetNp1()-1);
    }
    real_t gammaMax = sqrt(pMax*pMax + 1);
    real_t ximin = sqrt((gamma - 1) / (gamma + 1) * (gammaMax + 1) / (gammaMax - 1));

    if (abs(xi) < ximin)
        return 0;

    const real_t BA = this->grid->GetAvalancheCHBounceAverage(ir,i,j);
    return scaleFactor * preFactor * cutFactor * BA;
}

/**
 * Returns the source at grid point (ir,i,j).
 */
real_t AvalancheSourceCH::GetSourceFunction(len_t ir, len_t i, len_t j){
    real_t RESign;
    if (this->grid->GetMomentumGrid(ir)->GetP2(j) >= 0)
        RESign = 1;
    else
        RESign = -1;

    if(sourceMode == CH_SOURCE_MODE_FLUID){
        real_t S_f = 0.;
        real_t S, fre, fgrid;
        len_t offset = ir*n1[ir]*n2[ir];
        for (len_t ii=0; ii<n1[ir]; ii++){
            for (len_t jj=0; jj<n2[ir]; jj++){
                S = EvaluateCHSource(ir,ii,jj);
                len_t ind = offset + n1[ir]*jj + ii;
                
                real_t f_re_PAA = 0;
                
                int_t i_in = pIn_Indices[ind];
                if (htgridWithREgrid && i_in >= 0 && (len_t) i_in >= n1[ir]){
                    i_in -= n1[ir];
                    for(len_t j_int=0; j_int<runawayGrid->GetMomentumGrid(ir)->GetNp2(); j_int++){
                        if ((RESign >= 0 && this->runawayGrid->GetMomentumGrid(ir)->GetP2(j_int) >= 0)
                                || (RESign < 0 && this->runawayGrid->GetMomentumGrid(ir)->GetP2(j_int) < 0)){
                            fre = unknowns->GetUnknownData(id_fre)[nr*runawayGrid->GetMomentumGrid(ir)->GetNp1()*runawayGrid->GetMomentumGrid(ir)->GetNp2() + this->runawayGrid->GetMomentumGrid(ir)->GetNp1()*j_int + i_in];
                            f_re_PAA += fre * this->runawayGrid->GetMomentumGrid(ir)->GetDp2(j_int) / 2.;
                        }
                    }
                } else if (i_in >= 0) {
                    for(len_t j_int=0; j_int<n2[ir]; j_int++){
                        if ((RESign >= 0 && this->grid->GetMomentumGrid(ir)->GetP2(j_int) >= 0)
                                || (RESign < 0 && this->grid->GetMomentumGrid(ir)->GetP2(j_int) < 0)){
                            fgrid = unknowns->GetUnknownData(id_fgrid)[offset + n1[ir]*j_int + i_in];
                            f_re_PAA += fgrid * this->grid->GetMomentumGrid(ir)->GetDp2(j_int) / 2.;
                        }
                    }
                }
                S_f += S * f_re_PAA * this->grid->GetMomentumGrid(ir)->GetDp1(ii) * this->grid->GetMomentumGrid(ir)->GetDp2(jj);
            }
        }
        return S_f;
    }
    
    real_t S = EvaluateCHSource(ir,i,j);
    real_t fre, fgrid;
    len_t offset = ir*n1[ir]*n2[ir]; 
    len_t ind = offset + n1[ir]*j + i;
    
    real_t f_re_PAA = 0;
    
    int_t i_in = pIn_Indices[ind];
    if (htgridWithREgrid && i_in >= 0 && (len_t) i_in >= n1[ir]){
        i_in -= n1[ir];
        for(len_t j_int=0; j_int<runawayGrid->GetMomentumGrid(ir)->GetNp2(); j_int++){
            if ((RESign >= 0 && this->runawayGrid->GetMomentumGrid(ir)->GetP2(j_int) >= 0)
                    || (RESign < 0 && this->runawayGrid->GetMomentumGrid(ir)->GetP2(j_int) < 0)){
                fre = unknowns->GetUnknownData(id_fre)[ir*runawayGrid->GetMomentumGrid(ir)->GetNp1()*runawayGrid->GetMomentumGrid(ir)->GetNp2() + this->runawayGrid->GetMomentumGrid(ir)->GetNp1()*j_int + i_in];
                f_re_PAA += fre * this->runawayGrid->GetMomentumGrid(ir)->GetDp2(j_int) / 2.;
            }
        }
    } else if (i_in >= 0) {
        for(len_t j_int=0; j_int<n2[ir]; j_int++){
            if ((RESign >= 0 && this->grid->GetMomentumGrid(ir)->GetP2(j_int) >= 0)
                    || (RESign < 0 && this->grid->GetMomentumGrid(ir)->GetP2(j_int) < 0)){
                fgrid = unknowns->GetUnknownData(id_fgrid)[offset + n1[ir]*j_int + i_in];
                f_re_PAA += fgrid * this->grid->GetMomentumGrid(ir)->GetDp2(j_int) / 2.;
            }
        }
    }
    
    return S * f_re_PAA;
}

/**
 * Returns the source function at (ir,i,j) differentiated with respect to the unknown x_derivId at (ir,i,j)
 */
real_t AvalancheSourceCH::GetSourceFunctionJacobian(len_t ir, len_t i, len_t j, const len_t derivId){
    if (derivId==id_fgrid){
        real_t derivTerm = 0;
        len_t iter_j_gen = 0, j_gen;
        std::vector<len_t> f_xiIndices_i = f_xiIndices[ir*this->grid->GetMomentumGrid(ir)->GetNp1() + i];
        for (auto i_gen : f_pIndices[ir*this->grid->GetMomentumGrid(ir)->GetNp1() + i]) {
            j_gen = f_xiIndices_i[iter_j_gen];
            derivTerm += EvaluateCHSource(ir,i_gen,j_gen) * this->grid->GetMomentumGrid(ir)->GetDp2(j) / 2.;
            iter_j_gen++;
        }
        return derivTerm;
    }
    else if (htgridWithREgrid && derivId==id_fre) {
        real_t derivTerm = 0;
        len_t iter_j_gen = 0, j_gen;
        std::vector<len_t> f_xiIndices_i = f_re_xiIndices[ir*this->runawayGrid->GetMomentumGrid(ir)->GetNp1() + i];
        for (auto i_gen : f_re_pIndices[ir*this->runawayGrid->GetMomentumGrid(ir)->GetNp1() + i]) {
            j_gen = f_xiIndices_i[iter_j_gen];
            derivTerm += EvaluateCHSource(ir,i_gen,j_gen) * this->runawayGrid->GetMomentumGrid(ir)->GetDp2(j) / 2.;
            iter_j_gen++;
        }
        return derivTerm;
    }
    else
        return 0;
}

