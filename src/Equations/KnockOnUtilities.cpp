/**
 * Implementation of calculations that show up in kinetic models of
 * knock-on collisions of relativistic electrons (e.g. Møller scattering).
 *
 * The theory is described in detail in doc/notes/theory.tex,
 * and variable names throughout this file mirrors the corresponding
 * ones in the notes.
 *
 * There are a few recurring definitions throughout this implementation.
 * xi01 - is the incident electron pitch
 * [xi0_f1, xi0_f2] - is the knock-on pitch cell range over which we average the source.
 * xi_star - is the cosine of the true scattering angle (which is smeared out
 *   due to gyro-, orbit- and FVM averaging).
 */
#include "DREAM/Equations/KnockOnUtilities.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

#include "DREAM/DREAMException.hpp"

using namespace DREAM;

namespace {

constexpr real_t RELAXED_EPS = 5 * std::numeric_limits<real_t>::epsilon();


/*
 * RELAXED_EPS:
 * Used to guard against accumulated floating-point roundoff in chained
 * geometric and momentum-space transformations (e.g. xi0 → xi → xi1),
 * particularly near trapped–passing and endpoint boundaries.
 *
 * Chosen to be safely above machine epsilon while remaining negligible
 * on physical scales.
 *
 */

}  // anonymous namespace

/**
 * Analytic antiderivative of the gyro-averaged kinematic delta function.
 * Evaluates the variable "F" defined in the theory notes.
 */
real_t KnockOnUtilities::AnalyticDelta::EvaluateXiAntiderivative(
    real_t y, real_t xi1, real_t xi_star
) {
    const real_t argument=(y-xi1*xi_star)
        /sqrt((1-xi1)*(1+xi1)*(1-xi_star)*(1+xi_star));
    // Contacts sampled from either side can round the argument a few ulps
    // outside [-1,1]. Use the exact saturated antiderivative at those edges.
    return asin(std::clamp(argument,real_t(-1),real_t(1))) / M_PI;
}

/**
 * Analytic pitch integration of the gyro-averaged kinematic delta function,
 * from xi=t1 to xi=t2 (clipped to the kinematically available region [a-b, a+b],
 * in which case the values of F are limited to +/- 0.5).
 * Evaluates the variable "D" defined in the theory notes.
 */
real_t KnockOnUtilities::AnalyticDelta::EvaluateXiIntegral(
    real_t t1, real_t t2, real_t xi1, real_t xi_star
) {
    if (t2<=t1) return 0;
    real_t a = xi_star * xi1;
    real_t b2 = (1 - xi_star * xi_star) * (1 - xi1 * xi1);

    // handle singular (but analytically well-posed) case of
    // EvaluateXiAntiderivative when the denominator is 0
    if (b2 <= RELAXED_EPS) {
        real_t a = xi_star * xi1;
        return (t1 < a && a <= t2) ? 1.0 : 0.0;
    }

    real_t b = sqrt(b2);
    if (t2<=a-b || t1>=a+b) return 0;
    real_t F1 = -0.5;
    real_t F2 = 0.5;
    if (t2 < a + b) {
        F2 = EvaluateXiAntiderivative(t2, xi1, xi_star);
    }
    if (t1 > a - b) {
        F1 = EvaluateXiAntiderivative(t1, xi1, xi_star);
    }

    return F2 - F1;
}

/**
 * Compute the local bounds [t1, t2] obtained when changing from low-field side pitch [xi0_f1,
 * xi0_f2] in the FVM cell average to their local values xi(xi0, theta). The derivation of t1 and t2
 * is given in the theory notes.
 *
 * xi0Cutoff is the minimum low-field side pitch that can reach
 * this poloidal location (i.e. this value of BOverBmin).
 */
void KnockOnUtilities::AnalyticDelta::ComputeXiIntegrationBounds(
    real_t &t1, real_t &t2, real_t xi0_f1, real_t xi0_f2, real_t BOverBmin, real_t xi0Cutoff
) {
    // Clamp local pitch to zero where xi0 cannot reach this theta
    if (-RELAXED_EPS <= xi0_f1 && xi0_f1 <= xi0Cutoff) {
        t1 = 0;
    } else {
        t1 = FVM::MomentumGrid::XiFromXi0(xi0_f1, BOverBmin);
    }
    if (-xi0Cutoff <= xi0_f2 && xi0_f2 <= RELAXED_EPS) {
        t2 = 0;
    } else {
        t2 = FVM::MomentumGrid::XiFromXi0(xi0_f2, BOverBmin);
    }
}

/**
 * Evaluates the local limits [z1, z2] that xi1 spans when xi is swept between [t1, t2]
 * (e.g. when xi0 is cell-averaged over [xi0_f1, xi0_f2]), at fixed scattering angle cosine xi_star.
 * These limits are derived in the theory notes.
 */
void KnockOnUtilities::AnalyticDelta::ComputeXi1Bounds(
    real_t &z1, real_t &z2, real_t t1, real_t t2, real_t xi_star
) {
    if (t2 < -xi_star) {
        z1 = t2 * xi_star - sqrt(1 - t2 * t2) * sqrt(1 - xi_star * xi_star);
    } else if (t1 <= -xi_star && -xi_star <= t2) {
        z1 = -1;
    } else {
        z1 = t1 * xi_star - sqrt(1 - t1 * t1) * sqrt(1 - xi_star * xi_star);
    }

    if (t2 < xi_star) {
        z2 = t2 * xi_star + sqrt(1 - t2 * t2) * sqrt(1 - xi_star * xi_star);
    } else if (t1 <= xi_star && xi_star <= t2) {
        z2 = 1;
    } else {
        z2 = t1 * xi_star + sqrt(1 - t1 * t1) * sqrt(1 - xi_star * xi_star);
    }
}

/**
 * Evaluates the analytically gyro- and xi-cell-averaged kinematic
 * delta function at a fixed poloidal location, characterized
 * by its BOverBmin and xi0Cutoff values.
 *
 * In the theory notes, the returned value corresponds to the D_j function.
 *
 * Returns 0 if no contribution exists at this theta, and
 * handles the xi01 = +/-1 (CH/RP) limits carefully.
 */

real_t KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
    real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t BOverBmin, real_t xi_star, real_t xi0Cutoff,
    real_t &xi1_over_xi01
) {
    // perform the checks on whether the theta is viable
    // the "true" condition is z1 <= xi1 <= z2 for us to
    // take this theta into account. However, roundoff
    // errors when xi1 == z1, z2 can cause invalid values
    // to slip through and crash other parts of the calculation.
    // For the important edge case xi01 == +/- 1,
    // which define the Chiu-Harvey and Rosenbluth-Putvinski
    // approximations, we seek these bounding values, and therefore
    // employ an exacter constraint by inserting values explicitly.
    real_t t1;
    real_t t2;
    AnalyticDelta::ComputeXiIntegrationBounds(t1, t2, xi0_f1, xi0_f2, BOverBmin, xi0Cutoff);

    // handle the xi01 == 1 case (used in the approximate collision terms--RP and CH)
    if (fabs(fabs(xi01) - 1.0) < RELAXED_EPS) {
        xi1_over_xi01 = 1.0;
        // Use the same (lower,upper] ownership as the analytic gyro
        // integral, including the physical endpoint -1. Strict inequalities
        // on both faces would lose a pole source lying exactly on a face.
        const real_t outgoing = xi01 > 0 ? xi_star : -xi_star;
        if (!((t1 < outgoing && outgoing <= t2) || (outgoing == -1 && t1 == -1))) {
            return 0.0;
        }
        return 1.0 / (xi0_f2 - xi0_f1);
    }

    real_t z1;
    real_t z2;
    AnalyticDelta::ComputeXi1Bounds(z1, z2, t1, t2, xi_star);

    real_t xi1 = FVM::MomentumGrid::XiFromXi0(xi01, BOverBmin);

    if (fabs(xi01) > RELAXED_EPS) {
        xi1_over_xi01 = xi1 / xi01;
    } else {
        xi1_over_xi01 = 1.0;
    }
    if (!(z1 < xi1 && xi1 <= z2)) {
        return 0.0;
    }
    return AnalyticDelta::EvaluateXiIntegral(t1, t2, xi1, xi_star) / (xi0_f2 - xi0_f1);
}

/**
 * Evaluates the geometric quantities BOverBmin, Jacobian and xi0Cutoff at this poloidal location.
 * Returns true if this poloidal location can be reached both by some xi0 in [xi0_f1, xi0_f2] and by
 * xi01.
 */
bool KnockOnUtilities::CheckIfReachableAndSetGeometricQuantities(
    len_t ir, real_t theta, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t &BOverBmin,
    real_t &Jacobian, real_t &xi0Cutoff, FVM::FluxSurfaceAverager *fsa
) {
    real_t Bmin = fsa->GetBmin(ir, FVM::FLUXGRIDTYPE_DISTRIBUTION);
    real_t B_;
    real_t ROverR0_;
    real_t NablaR2_;
    fsa->GeometricQuantitiesAtTheta(
        ir, theta, B_, Jacobian, ROverR0_, NablaR2_, FVM::FLUXGRIDTYPE_DISTRIBUTION
    );
    BOverBmin = Bmin == 0 && B_ == 0 ? 1 : B_ / Bmin;
    // Clamp BOverBmin to avoid issues with numerical geometry
    BOverBmin = std::max(BOverBmin, 1.0);

    xi0Cutoff = sqrt(1 - 1 / BOverBmin);

    if (xi0Cutoff >= fabs(xi01) && BOverBmin > 1) {
        // xi01 orbit cannot reach this poloidal location => return 0.
        return false;
    }
    if (xi0Cutoff >= std::max(fabs(xi0_f1), fabs(xi0_f2))) {
        // no xi0 orbit in [xi0_f1, xi0_f2] can reach this poloidal location => return 0.
        return false;
    }
    return true;
}

/**
 * Orbit-dependent integrand of the analytically gyro- and FVM-averaged
 * kinematic delta function.
 *
 * Returns zero whenever kinematic or trapping constraints are violated
 * at the given poloidal angle theta.
 */
real_t KnockOnUtilities::OrbitDeltaIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::RadialGrid *rg
) {
    real_t BOverBmin;
    real_t Jacobian;
    real_t xi0Cutoff;
    bool isReachable = CheckIfReachableAndSetGeometricQuantities(
        ir, theta, xi0_f1, xi0_f2, xi01, BOverBmin, Jacobian, xi0Cutoff,
        rg->GetFluxSurfaceAverager()
    );
    if (!isReachable) {
        return 0.0;
    }

    real_t xi1_over_xi01;
    real_t D = AnalyticDelta::EvaluateLocalContribution(
        xi0_f1, xi0_f2, xi01, BOverBmin, xi_star, xi0Cutoff, xi1_over_xi01
    );

    real_t sqrtG = Jacobian * FVM::MomentumGrid::evaluatePXiMetricOverP2(xi1_over_xi01, BOverBmin);

    return sqrtG * D;
}

real_t KnockOnUtilities::OrbitNumeratorIntegrand(
    real_t theta, len_t ir, real_t xi0_f1, real_t xi0_f2, real_t xi01, real_t xi_star,
    const FVM::RadialGrid *rg
) {
    // OrbitDeltaIntegrand contains the metric over p1^2 and the outgoing
    // pitch-cell average. Supply the remaining toroidal factor, without ever
    // dividing by the (singular at the separatrix) pointwise primary volume.
    return 2*M_PI*OrbitDeltaIntegrand(theta,ir,xi0_f1,xi0_f2,xi01,xi_star,rg);
}

namespace {
namespace Q = KnockOnUtilities::Quadrature;

struct OrbitPiece {
    Q::Interval cell;
    len_t row;
    real_t weight;
    bool mirror;
};

// Split before folding so that each piece has one trapped multiplicity.
// Weights convert a subcell average back to its original row average.
void AddTrappedPieces(
    std::vector<OrbitPiece>& pieces, real_t lo, real_t hi, len_t row,
    real_t u, real_t xiT
) {
    const auto faces=Q::SplitPoints({lo,hi},{-xiT,0,xiT});
    const bool primaryTrapped=u>=0 && u<xiT;
    for (len_t i=1; i<faces.size(); ++i) {
        const real_t a=faces[i-1], b=faces[i];
        if (a>=-xiT && b<=0) continue;
        const bool secondaryTrapped=a>=0 && b<=xiT;
        pieces.push_back({{a,b},row,
            (b-a)/(hi-lo)*(primaryTrapped && secondaryTrapped ? 2 : 1),
            primaryTrapped || secondaryTrapped});
    }
}

// Both trapped signs use the same contacts, metric and geometric quantities.
real_t CombinedOrbitIntegrand(
    real_t theta, len_t ir, Q::Interval cell, real_t u, real_t chi,
    bool mirror, const FVM::RadialGrid *rg
) {
    real_t B, J, cutoff;
    if (!KnockOnUtilities::CheckIfReachableAndSetGeometricQuantities(
        ir,theta,cell.lower,cell.upper,u,B,J,cutoff,rg->GetFluxSurfaceAverager()
    )) return 0;
    real_t ratio;
    real_t D=KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
        cell.lower,cell.upper,u,B,chi,cutoff,ratio
    );
    if (mirror) {
        real_t unused;
        D+=KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
            cell.lower,cell.upper,u,B,-chi,cutoff,unused
        );
    }
    return 2*M_PI*J*FVM::MomentumGrid::evaluatePXiMetricOverP2(ratio,B)*D;
}

void IntegratePieces(
    len_t ir, real_t chi, real_t u, const std::vector<OrbitPiece>& pieces,
    const FVM::RadialGrid *rg, Q::OrbitGeometry& geometry,
    const Q::GaussLegendreRule& rule, real_t *out
) {
    if (pieces.empty()) return;
    if (geometry.IsConstantField()) {
        const real_t h=2*M_PI*rg->GetVpVol(ir);
        for (const auto& piece : pieces) {
            real_t ratio;
            real_t D=KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
                piece.cell.lower,piece.cell.upper,u,1,chi,0,ratio
            );
            if (piece.mirror)
                D+=KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
                    piece.cell.lower,piece.cell.upper,u,1,-chi,0,ratio
                );
            out[piece.row]+=piece.weight*h*D;
        }
        return;
    }
    std::vector<Q::Interval> cells;
    cells.reserve(pieces.size());
    for (const auto& piece : pieces) cells.push_back(piece.cell);
    const auto panels=Q::PoloidalSplitPoints(geometry,u,chi,cells);
    for (len_t j=0; j<pieces.size(); ++j) {
        const auto& piece=pieces[j];
        real_t integral=0;
        for (len_t k=1; k<panels[j].size(); ++k)
            for (const auto& point : rule.MappedPoints({panels[j][k-1],panels[j][k]}))
                integral+=point.weight*CombinedOrbitIntegrand(
                    point.x,ir,piece.cell,u,chi,piece.mirror,rg
                );
        if (!std::isfinite(integral) || integral<0)
            throw DREAMException(
                "Knock-on orbit quadrature produced numerator %.17g at ir=%ld, "
                "u=%.17g, chi=%.17g, cell=[%.17g,%.17g].",
                integral,(long)ir,u,chi,piece.cell.lower,piece.cell.upper
            );
        out[piece.row]+=piece.weight*integral;
    }
}

// A target piece carries a density per ORIGINAL outgoing row width. This
// includes a void boundary cell's passing sliver reassigned to a different row.
struct TargetPiece {
    Q::Interval cell;
    len_t row;
    real_t weight;
    bool trapped;
};
void AddTargetPieces(
    std::vector<TargetPiece>& pieces, Q::Interval cell, len_t row,
    real_t rowWidth, real_t xiT
) {
    const auto cuts=Q::SplitPoints(cell,{-xiT,0,xiT});
    for (len_t k=1; k<cuts.size(); ++k) {
        const real_t a=cuts[k-1], b=cuts[k];
        if (a>=-xiT && b<=0) continue;
        pieces.push_back({{a,b},row,(b-a)/rowWidth,a>=0 && b<=xiT});
    }
}

void IntegratePrimaryCell(
    len_t ir, real_t chi, Q::Interval primary, const std::vector<TargetPiece>& targets,
    const FVM::RadialGrid *rg, Q::OrbitGeometry& geometry,
    const Q::GaussLegendreRule& poloidalRule, const Q::GaussLegendreRule& pitchRule,
    real_t *out
) {
    const real_t betaMax=geometry.GetBetaMax(), T=std::sqrt(betaMax);
    // Each target piece gets its own outer partition. Batch coincident nodes
    // locally so repeated field levels are still inverted once. No cache
    // survives this primary column, and no row pays for another row's cuts.
    std::map<real_t,std::vector<OrbitPiece>> samples;
    for (const auto& target : targets) {
        const auto cuts=Q::PrimaryPitchSplitPoints(
            primary,chi,betaMax,{target.cell.lower,target.cell.upper}
        );
        for (len_t k=1; k<cuts.size(); ++k) {
            if (cuts[k-1]>=-T && cuts[k]<=0) continue;
            for (const auto& point : pitchRule.MappedPoints({cuts[k-1],cuts[k]})) {
                const real_t u=point.x;
                // A sub-ulp panel can round onto the separatrix itself, whose
                // pointwise measure is infinite. Such a sample is unresolved
                // by this field representation; never evaluate that infinity.
                if (betaMax>0 && (std::abs(u)==T || u*u==betaMax)) continue;
                const bool trapped=u>=0 && u<T;
                samples[u].push_back({target.cell,target.row,
                    point.weight*target.weight*(trapped && target.trapped ? 2 : 1),
                    trapped || target.trapped});
            }
        }
    }
    for (const auto& sample : samples)
        IntegratePieces(ir,chi,sample.first,sample.second,rg,geometry,poloidalRule,out);
}

// Raw cell numerator divided by DREAM's discrete primary measure Hgrid.
// It is NOT divided by the pointwise orbit volume at any quadrature node.
void EvaluateGridRows(
    len_t ir, real_t chi, len_t l, const FVM::Grid *gridK, const FVM::Grid *gridP,
    len_t first, len_t last, Q::OrbitGeometry& geometry,
    const Q::GaussLegendreRule& poloidalRule, const Q::GaussLegendreRule& pitchRule,
    real_t *out
) {
    std::fill(out,out+last-first,0);
    const auto *mgP=gridP->GetMomentumGrid(ir), *mgK=gridK->GetMomentumGrid(ir);
    const real_t Vp1=gridP->GetVpOverP2AtZero(ir)[l];
    // Native void columns stay void, including a negative-passing sliver
    // belonging to the discarded cell crossing -xiT.
    if (Vp1==0) return;
    const real_t Hgrid=mgP->GetDp2(l)*Vp1;
    if (!(Hgrid>0) || !std::isfinite(Hgrid))
        throw DREAMException("Knock-on primary cell has invalid volume.");
    const real_t T=std::sqrt(geometry.GetBetaMax());
    std::vector<TargetPiece> targets;
    for (len_t j=first; j<last; ++j) {
        if (gridK->GetVpOverP2AtZero(ir)[j]==0) continue;
        const real_t lo=mgK->GetP2_f(j), hi=mgK->GetP2_f(j+1), width=mgK->GetDp2(j);
        AddTargetPieces(targets,{lo,hi},j-first,width,T);
        if (T>0 && lo<=T && T<hi) {
            // Restore only the passing part of a native void boundary row.
            // Divide by the RECEIVING row width; this is essential on unlike
            // positive/negative pitch grids. Never transfer primary inventory.
            for (len_t k=0; k<mgK->GetNp2(); ++k)
                if (mgK->GetP2_f(k)<-T && mgK->GetP2_f(k+1)>-T
                    && gridK->GetVpOverP2AtZero(ir)[k]==0)
                    AddTargetPieces(targets,{mgK->GetP2_f(k),-T},j-first,width,T);
        }
    }
    IntegratePrimaryCell(ir,chi,{mgP->GetP2_f(l),mgP->GetP2_f(l+1)},targets,
        gridP->GetRadialGrid(),geometry,poloidalRule,pitchRule,out);
    for (len_t j=0; j<last-first; ++j) out[j]/=Hgrid;
}
} // namespace

// Explicit-pitch, unnormalized interface: no grid-centre substitution or
// primary-volume division. Also used by the independent numerical tests.
real_t KnockOnUtilities::EvaluateOrbitNumerator(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry& geometry,
    const Quadrature::GaussLegendreRule& rule, bool mirror
) {
    if (!(xi0_f2>xi0_f1))
        throw DREAMException("Knock-on orbit quadrature needs a nonempty outgoing cell.");
    real_t result=0;
    IntegratePieces(ir,xi_star,xi01,{{{xi0_f1,xi0_f2},0,1,mirror}},rg,geometry,rule,&result);
    return result;
}

real_t KnockOnUtilities::EvaluatePitchCellNumerator(
    len_t ir, real_t xi_star, Quadrature::Interval primary, Quadrature::Interval secondary,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry& geometry,
    const Quadrature::GaussLegendreRule& poloidalRule,
    const Quadrature::GaussLegendreRule& pitchRule
) {
    if (!(primary.upper>primary.lower) || !(secondary.upper>secondary.lower))
        throw DREAMException("Knock-on cell integral requires nonempty pitch intervals.");
    std::vector<TargetPiece> targets;
    AddTargetPieces(targets,secondary,0,secondary.upper-secondary.lower,std::sqrt(geometry.GetBetaMax()));
    real_t result=0;
    IntegratePrimaryCell(ir,xi_star,primary,targets,rg,geometry,poloidalRule,pitchRule,&result);
    return result;
}

real_t KnockOnUtilities::EvaluateOrbitAveragedDelta(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1,
    const FVM::RadialGrid *rg, len_t poloidalOrder
) {
    Quadrature::OrbitGeometry geometry(rg->GetFluxSurfaceAverager(),ir);
    Quadrature::GaussLegendreRule rule(poloidalOrder);
    if (!(Vp1>0)) throw DREAMException("Knock-on orbit average requires positive primary volume.");
    return EvaluateOrbitNumerator(ir,xi_star,xi01,xi0_f1,xi0_f2,rg,geometry,rule)/Vp1;
}

real_t KnockOnUtilities::EvaluateOrbitAveragedDeltaWithTrappingCorrection(
    len_t ir, real_t xi_star, real_t xi01, real_t xi0_f1, real_t xi0_f2, real_t Vp1,
    const FVM::RadialGrid *rg, len_t poloidalOrder
) {
    Quadrature::OrbitGeometry geometry(rg->GetFluxSurfaceAverager(),ir);
    Quadrature::GaussLegendreRule rule(poloidalOrder);
    if (!(Vp1>0) || !(xi0_f2>xi0_f1))
        throw DREAMException("Knock-on orbit average requires positive volume and cell width.");
    std::vector<OrbitPiece> pieces;
    AddTrappedPieces(pieces,xi0_f1,xi0_f2,0,xi01,rg->GetXi0TrappedBoundary(ir));
    real_t result=0;
    IntegratePieces(ir,xi_star,xi01,pieces,rg,geometry,rule,&result);
    return result/Vp1;
}

real_t KnockOnUtilities::EvaluateDeltaMatrixElementOnGrid(
    len_t ir, real_t xi_star, len_t j, len_t l, const FVM::Grid *grid_knockon,
    const FVM::Grid *grid_primary, len_t poloidalOrder, len_t pitchOrder
) {
    Quadrature::OrbitGeometry geometry(grid_primary->GetRadialGrid()->GetFluxSurfaceAverager(),ir);
    Quadrature::GaussLegendreRule rule(poloidalOrder), pitchRule(pitchOrder);
    real_t result;
    EvaluateGridRows(ir,xi_star,l,grid_knockon,grid_primary,j,j+1,geometry,rule,pitchRule,&result);
    return result;
}

real_t KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, len_t poloidalOrder, len_t pitchOrder
) {
    Quadrature::OrbitGeometry geometry(grid_primary->GetRadialGrid()->GetFluxSurfaceAverager(),ir);
    Quadrature::GaussLegendreRule rule(poloidalOrder), pitchRule(pitchOrder);
    return SetDeltaMatrixColumnOnGrid(ir,xi_star,l,grid_knockon,grid_primary,deltaCol,geometry,rule,pitchRule);
}

real_t KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, Quadrature::OrbitGeometry& geometry, const Quadrature::GaussLegendreRule& rule,
    const Quadrature::GaussLegendreRule& pitchRule
) {
    const auto *mg=grid_knockon->GetMomentumGrid(ir);
    const len_t n=mg->GetNp2();
    EvaluateGridRows(ir,xi_star,l,grid_knockon,grid_primary,0,n,geometry,rule,pitchRule,deltaCol);
    if (grid_primary->GetVpOverP2AtZero(ir)[l]==0) return 0;
    real_t norm=0;
    for (len_t j=0; j<n; ++j) norm+=mg->GetDp2(j)*deltaCol[j];
    if (!(norm>0) || !std::isfinite(norm))
        throw DREAMException(
            "Knock-on quadrature expected a positive column at ir=%ld, l=%ld, xi_star=%.16g.",
            (long)ir,(long)l,xi_star
        );
    for (len_t j=0; j<n; ++j) deltaCol[j]/=norm;
    return norm;
}

/**
 * Evaluates the momentum-space knock-on kernel S_{ik} on a DREAM grid.
 *
 * Assumes a p-xi grid that is identical at all radii.
 *
 * Uses p1 on primary grid centers at index k and p_{i-1/2}, p_{i+1/2}
 * on the knock-on momentum grid at index i.
 *
 * This function is a grid-aware convenience wrapper around
 * Kinematics::EvaluateMollerFlux(), including the primary speed factor v1,
 * with kinematic clamping to the relativistically allowed knock-on domain.
 *
 * Corresponds to the matrix element S_{ik} in the theory notes.
 */
real_t KnockOnUtilities::EvaluateMollerFluxMatrixElementOnGrid(
    len_t i, len_t k, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t pCutoff
) {
    if (!(pCutoff > RELAXED_EPS)) {
        throw DREAMException(
            "EvaluateMollerFluxMatrixElementOnGrid(): invalid pCutoff=%.16g (must satisfy pCutoff "
            "> %.3g). "
            "This is required because the Møller flux antiderivative is singular at p=0.",
            pCutoff, RELAXED_EPS
        );
    }

    real_t p_f1 = grid_knockon->GetMomentumGrid(0)->GetP1_f(i);
    real_t p_f2 = grid_knockon->GetMomentumGrid(0)->GetP1_f(i + 1);

    if (p_f2 <= pCutoff) {
        // entire momentum cell below the admissible interval
        return 0.0;
    }

    real_t p1 = grid_primary->GetMomentumGrid(0)->GetP1(k);
    real_t pMax = KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1);

    if (p_f1 > pMax) {
        // entire momentum cell above the admissible interval
        return 0.0;
    }

    // Clamp range to [pCutoff, pMax]
    real_t pLower = std::max(p_f1, pCutoff);
    real_t pUpper = std::min(p_f2, pMax);
    if (!(pUpper > pLower + RELAXED_EPS)) {
        return 0.0;
    }

    real_t S2 = KnockOnUtilities::Kinematics::EvaluateMollerFlux(pUpper, p1);
    real_t S1 = KnockOnUtilities::Kinematics::EvaluateMollerFlux(pLower, p1);
    real_t v1 = p1 / sqrt(1 + p1 * p1);
    return v1 * (S2 - S1) / (p_f2 - p_f1);
}

/**
 * Evaluate the total Møller flux on the entire knock-on grid
 * This is an O(1) shortcut to evaluate the commonly appearing
 *   sum_i dp_i * EvaluateMollerFluxMatrixElementOnGrid(i,k,...)
 * that appears when the total production rate of knock-ons due
 * to a primary distribution is evaluated.
 */
real_t KnockOnUtilities::EvaluateMollerFluxIntegratedOverKnockonGrid(
    len_t k, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary, real_t pCutoff
) {
    if (!(pCutoff > RELAXED_EPS)) {
        throw DREAMException(
            "EvaluateMollerFluxIntegratedOverKnockonGrid(): invalid pCutoff=%.16g (must satisfy "
            "pCutoff > %.3g).",
            pCutoff, RELAXED_EPS
        );
    }

    const auto *mgK = grid_knockon->GetMomentumGrid(0);
    const auto *mgP = grid_primary->GetMomentumGrid(0);

    const real_t p1 = mgP->GetP1(k);
    const real_t v1 = p1 / sqrt(1 + p1 * p1);

    const real_t pMax = KnockOnUtilities::Kinematics::MaximumKnockOnMomentum(p1);

    // Knock-on grid extent
    const real_t pBottom = mgK->GetP1_f(0);
    const real_t pTop = mgK->GetP1_f(mgK->GetNp1());

    const real_t plo = std::max(pCutoff, pBottom);
    const real_t phi = std::min(pMax, pTop);

    if (!(phi > plo + RELAXED_EPS)) {
        return 0.0;
    }

    const real_t Fhi = KnockOnUtilities::Kinematics::EvaluateMollerFlux(phi, p1);
    const real_t Flo = KnockOnUtilities::Kinematics::EvaluateMollerFlux(plo, p1);

    // This equals sum_i dp_i * EvaluateMollerFluxMatrixElementOnGrid(i,k,...)
    return v1 * (Fhi - Flo);
}
