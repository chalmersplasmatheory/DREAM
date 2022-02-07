#include "DREAM/Equations/Kinetic/BraamsKarneyRHS.hpp"

using namespace DREAM;

static bool DoIntegralBC(len_t i, len_t j, len_t np1, len_t np2) {
    (void)j, (void)np2;
    return i == np1 - 1;
}

BraamsKarneyRHS::BraamsKarneyRHS(FVM::Grid *grid)
    : EquationTerm(grid) {
    GridRebuilt();
}

len_t BraamsKarneyRHS::GetNumberOfNonZerosPerRow() const {
    return 1;
}

bool BraamsKarneyRHS::SetJacobianBlock(const len_t unknId, const len_t derivId, FVM::Matrix *jac, const real_t *) {
    if (derivId == unknId) {
        this->SetMatrixElements(jac, nullptr);
        return true;
    }

    return false;
}

void BraamsKarneyRHS::SetMatrixElements(FVM::Matrix *mat, real_t*) {
    SetElements([&](auto idx, auto idx_p, auto V) {
                    mat->SetElement(idx, idx_p, V);
                });
}

void BraamsKarneyRHS::SetVectorElements(real_t *vec, const real_t *f) {
    SetElements([&](auto idx, auto idx_p, auto V) {
                    vec[idx] += f[idx_p] * V;
                });
}

template<typename F1>
void BraamsKarneyRHS::SetElements(F1 X) {
    const len_t nr  = grid->GetNr();

    len_t offset = 0;
    for (len_t ir = 0; ir < nr; ir++) {
        const FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        const len_t np1 = mg->GetNp1();
        const len_t np2 = mg->GetNp2();

        for (len_t i = 0; i < np1; i++) {
            for (len_t j = 0; j < np2; j++) {
                len_t idx = offset + j * np1 + i;

                if (!DoIntegralBC(i, j, np1, np2)) {
                    if (i == 0 || j == 0 || j == np2 - 1) {
                        // RHS should be zero, since these
                        // represent d Psi / d p = 0 and d Psi / d xi = 0.
                    } else {
                        X(idx, idx, -1);
                    }
                }
            }
        }

        offset += np1 * np2;
    }
}
