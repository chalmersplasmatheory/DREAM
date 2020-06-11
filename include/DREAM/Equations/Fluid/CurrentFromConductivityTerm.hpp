#include "FVM/Equation/DiagonalComplexTerm.hpp"
#include "DREAM/Equations/RunawayFluid.hpp"
#include "DREAM/IonHandler.hpp"
/**
 * Implementation of a class which represents the sigma*E contribution to the ohmic current equation.
 */
namespace DREAM {
    class CurrentFromConductivityTerm : public FVM::DiagonalComplexTerm {
    private:
        RunawayFluid *REFluid;
        IonHandler *ionHandler;
    protected:
        virtual bool TermDependsOnUnknowns() override {return true;}
        virtual void SetDiffWeights(len_t derivId, len_t nMultiples) override {
            
            real_t *dSigma = REFluid->evaluatePartialContributionConductivity(ionHandler->evaluateZeff(),derivId);

            len_t offset = 0;
            for(len_t n = 0; n<nMultiples; n++){
                for (len_t ir = 0; ir < nr; ir++){
                    //real_t w=0;
                    real_t w0 = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir));
                    for(len_t i = 0; i < n1[ir]; i++)
                        for(len_t j = 0; j < n2[ir]; j++)
                            diffWeights[offset + n1[ir]*j + i] = w0*dSigma[offset + n1[ir]*j + i];
                    offset += n1[ir]*n2[ir];
                }
            }
        }
    public:
        CurrentFromConductivityTerm(FVM::Grid* g, FVM::UnknownQuantityHandler *u, RunawayFluid *ref, IonHandler *ih) 
            : FVM::DiagonalComplexTerm(g,u), REFluid(ref), ionHandler(ih)
        {
            AddUnknownForJacobian(unknowns,unknowns->GetUnknownID(OptionConstants::UQTY_T_COLD));
        }

        virtual void SetWeights() override {
            len_t offset = 0;
            for (len_t ir = 0; ir < nr; ir++){
                //real_t w=0;
                real_t w = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir)) * REFluid->evaluateElectricalConductivity(ir, ionHandler->evaluateZeff(ir));
                for(len_t i = 0; i < n1[ir]; i++)
                    for(len_t j = 0; j < n2[ir]; j++)
                        weights[offset + n1[ir]*j + i] = w;
                offset += n1[ir]*n2[ir];
            }
        }
    };
}