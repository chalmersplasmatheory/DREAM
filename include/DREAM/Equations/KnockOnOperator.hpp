#ifndef _DREAM_EQUATIONS_KNOCK_ON_OPERATOR_HPP
#define _DREAM_EQUATIONS_KNOCK_ON_OPERATOR_HPP

namespace DREAM {
class KnockOnOperator;
}

#include "DREAM/DREAMException.hpp"
#include "FVM/Grid/Grid.hpp"
#include <cmath>

namespace DREAM {
class KnockOnOperator {
  private:
    const FVM::Grid *grid;
    real_t p_cutoff;
    len_t n_xi_stars_tabulate;

    // tabulated values of size Nr x n_xi_stars_tabulate x (Nxi*Nxi)
    real_t ***deltaTable = nullptr;
    real_t *xi_stars_tab = nullptr;
    void Allocate();
    void Deallocate();

  public:
    KnockOnOperator(
        const FVM::Grid *grid, real_t p_cutoff, len_t n_xi_stars_tabulate = 50,
        len_t n_points_integral = 80
    );
    ~KnockOnOperator();

    bool GridRebuilt();

    real_t EvaluateDelta(len_t ir, real_t xi_star, len_t j, len_t l);

    void TabulateDeltaOnGrid();

    len_t n_points_integral;
};
} // namespace DREAM

#endif /*_DREAM_EQUATIONS_KNOCK_ON_OPERATOR_HPP*/
