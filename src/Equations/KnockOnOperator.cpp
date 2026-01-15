/**
 * Implements the knock-on collision operator for binary large-angle collisions.
 */
#include "DREAM/Equations/KnockOnOperator.hpp"
#include "DREAM/Equations/KnockOnUtilities.hpp"
#include <cmath>
#include <limits>

using namespace DREAM;

KnockOnOperator::KnockOnOperator(
    const FVM::Grid *grid, real_t p_cutoff, len_t n_xi_stars_tabulate, len_t n_points_integral
)
    : grid(grid), p_cutoff(p_cutoff), n_xi_stars_tabulate(n_xi_stars_tabulate),
      n_points_integral(n_points_integral) {
    Allocate();
}

KnockOnOperator::~KnockOnOperator() { Deallocate(); }

void KnockOnOperator::Deallocate() {
    len_t Nr = grid->GetNr();
    for (len_t ir = 0; ir < Nr; ir++) {
        // XXX: assume p-xi grid
        deltaTable[ir] = new real_t *[n_xi_stars_tabulate];
        for (len_t m = 0; m < n_xi_stars_tabulate; m++) {
            delete[] deltaTable[ir][m];
        }
        delete[] deltaTable[ir];
    }
    delete[] deltaTable;
    delete[] xi_stars_tab;
}

void KnockOnOperator::Allocate() {
    len_t Nr = grid->GetNr();
    deltaTable = new real_t **[Nr];
    for (len_t ir = 0; ir < Nr; ir++) {
        len_t Nxi = grid->GetNp2(ir);
        deltaTable[ir] = new real_t *[n_xi_stars_tabulate];
        for (len_t m = 0; m < n_xi_stars_tabulate; m++) {
            deltaTable[ir][m] = new real_t[Nxi * Nxi];
        }
    }

    // XXX: assume same momentum grid at all radii, and p-xi grid
    real_t p1_max = grid->GetMomentumGrid(0)->GetP(grid->GetNp1(0) - 1, grid->GetNp2(0) - 1);
    real_t g1_max = sqrt(1 + p1_max * p1_max);
    real_t g_max = (g1_max + 1) / 2;
    real_t p_max = sqrt(g_max * g_max - 1);

    real_t xi_star_min = KnockOnUtilities::evaluateXiStar(p_cutoff, p1_max);
    real_t xi_star_max = KnockOnUtilities::evaluateXiStar(p_max, p1_max);
    real_t dxi = (xi_star_max - xi_star_min) / (n_xi_stars_tabulate - 1);
    xi_stars_tab = new real_t[n_xi_stars_tabulate];
    for (len_t i = 0; i < n_xi_stars_tabulate; i++) {
        xi_stars_tab[i] = xi_star_min + i * dxi;
    }
}

bool KnockOnOperator::GridRebuilt() {
    Deallocate();
    Allocate();
    return true;
}

/**
 * Evaluate delta on the grid, assuming a p-xi grid.
 * Evaluates xi01 on grid centers at index l and xi0_{j-1/2}, xi0_{j+1/2} at index j.
 */
real_t KnockOnOperator::EvaluateDelta(len_t ir, real_t xi_star, len_t j, len_t l) {
    // no contributions to the negative trapped region -
    // instead (below) these contributions are added to the positive trapped

    FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
    real_t xi01 = mg->GetXi0(0, l);
    real_t xi0_f1 = mg->GetXi0_f2(0, j);
    real_t xi0_f2 = mg->GetXi0_f2(0, j + 1);
    real_t Vp1 = grid->GetVpOverP2AtZero(ir)[l];

    // no contribution for unphysical parts of the xi01 phase space
    if (!Vp1)
        return 0;

    // Set the outer theta interval as the smallest of the orbits reached by xi_{01} and
    // [xi_{0,j-1/2}, xi_{0, j+1/2}]
    real_t theta1, theta2;
    KnockOnUtilities::estimateBoundingTheta(ir, j, l, theta1, theta2, grid);

    real_t delta = KnockOnUtilities::EvaluateDelta(
        ir, xi_star, xi01, xi0_f1, xi0_f2, Vp1, theta1, theta2, n_points_integral, grid
    );
    return delta;
}

void KnockOnOperator::TabulateDeltaOnGrid() {
    len_t Nr = grid->GetNr();
    for (len_t ir = 0; ir < Nr; ir++) {
        len_t Nxi = grid->GetNp2(ir);
        FVM::MomentumGrid *mg = grid->GetMomentumGrid(ir);
        real_t *deltaRow = new real_t[Nxi];
        for (len_t m = 0; m < n_xi_stars_tabulate; m++) {
            for (len_t l = 0; l < Nxi; l++) {
                // evaluate a full j-row, and then renormalize
                // the result to the known exact relation:
                //     \sum_j (xi_{j+1/2} - xi_{j-1/2}) * delta_{jl} = 1
                // to compensate for (small) quadrature errors.
                real_t deltaSum = 0;
                for (len_t j = 0; j < Nxi; j++) {
                    real_t delta = EvaluateDelta(ir, xi_stars_tab[m], j, l);
                    deltaRow[j] = delta;
                    real_t dxi = mg->GetDp2(j);
                    deltaSum += dxi * delta;
                }
                // set normalized values in the table
                for (len_t j = 0; j < Nxi; j++) {
                    len_t idx = l * Nxi + j;
                    deltaTable[ir][m][idx] = deltaRow[j] / deltaSum;
                }
            }
        }
        delete[] deltaRow;
    }
}
