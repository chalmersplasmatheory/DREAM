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

namespace Quadrature = KnockOnUtilities::Quadrature;

// Outgoing subcell to integrate at one primary-pitch sample. The weight
// converts its contribution to the stored row and includes multiplicity.
struct OrbitPiece {
    Quadrature::Interval cell;
    len_t row;
    real_t weight;
    bool mirror;
};

// Split a fixed-pitch orbit integral before applying trapped multiplicities.
void AddTrappedPieces(
    std::vector<OrbitPiece> &pieces, real_t lo, real_t hi, len_t row,
    real_t primaryPitch, real_t xiT
) {
    const auto faces = Quadrature::SplitPoints({lo, hi}, {-xiT, 0, xiT});
    const bool primaryTrapped = primaryPitch >= 0 && primaryPitch < xiT;
    for (len_t i = 1; i < faces.size(); ++i) {
        const real_t lower = faces[i - 1];
        const real_t upper = faces[i];
        if (lower >= -xiT && upper <= 0)
            continue;

        const bool secondaryTrapped = lower >= 0 && upper <= xiT;
        const int multiplicity = primaryTrapped && secondaryTrapped ? 2 : 1;
        pieces.push_back({
            {lower, upper}, row, (upper - lower) / (hi - lo) * multiplicity,
            primaryTrapped || secondaryTrapped
        });
    }
}

// Both trapped signs use the same contacts, metric and geometric quantities.
real_t CombinedOrbitIntegrand(
    real_t theta, len_t ir, Quadrature::Interval cell, real_t primaryPitch, real_t xiStar,
    bool mirror, const FVM::RadialGrid *rg
) {
    real_t BOverBmin;
    real_t jacobian;
    real_t xi0Cutoff;
    const bool reachable = KnockOnUtilities::CheckIfReachableAndSetGeometricQuantities(
        ir, theta, cell.lower, cell.upper, primaryPitch, BOverBmin, jacobian, xi0Cutoff,
        rg->GetFluxSurfaceAverager()
    );
    if (!reachable)
        return 0;

    real_t xi1OverXi01;
    real_t delta = KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
        cell.lower, cell.upper, primaryPitch, BOverBmin, xiStar, xi0Cutoff, xi1OverXi01
    );

    if (mirror) {
        real_t unused;
        delta += KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
            cell.lower, cell.upper, primaryPitch, BOverBmin, -xiStar, xi0Cutoff, unused
        );
    }

    return 2 * M_PI * jacobian
        * FVM::MomentumGrid::evaluatePXiMetricOverP2(xi1OverXi01, BOverBmin) * delta;
}

/**
 * Accumulate the poloidal numerators for one primary pitch. Each outgoing
 * piece has its own panels, with common field levels inverted as a batch.
 */
void IntegratePieces(
    len_t ir, real_t xiStar, real_t primaryPitch, const std::vector<OrbitPiece> &pieces,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry &geometry,
    const Quadrature::GaussLegendreRule &poloidalRule, real_t *numerators
) {
    if (pieces.empty())
        return;

    if (geometry.IsConstantField()) {
        const real_t orbitMeasure = 2 * M_PI * rg->GetVpVol(ir);
        for (const auto &piece : pieces) {
            real_t ratio;
            real_t delta = KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
                piece.cell.lower, piece.cell.upper, primaryPitch, 1, xiStar, 0, ratio
            );
            if (piece.mirror) {
                delta += KnockOnUtilities::AnalyticDelta::EvaluateLocalContribution(
                    piece.cell.lower, piece.cell.upper, primaryPitch, 1, -xiStar, 0, ratio
                );
            }
            numerators[piece.row] += piece.weight * orbitMeasure * delta;
        }
        return;
    }

    std::vector<Quadrature::Interval> cells;
    cells.reserve(pieces.size());
    for (const auto &piece : pieces) {
        cells.push_back(piece.cell);
    }

    const auto panels = Quadrature::PoloidalSplitPoints(geometry, primaryPitch, xiStar, cells);
    for (len_t j = 0; j < pieces.size(); ++j) {
        const auto &piece = pieces[j];
        const auto &cuts = panels[j];
        real_t integral = 0;
        for (len_t k = 1; k < cuts.size(); ++k) {
            const Quadrature::Interval panel{cuts[k - 1], cuts[k]};
            for (const auto &point : poloidalRule.MappedPoints(panel)) {
                integral += point.weight * CombinedOrbitIntegrand(
                    point.x, ir, piece.cell, primaryPitch, xiStar, piece.mirror, rg
                );
            }
        }
        if (!std::isfinite(integral) || integral < 0) {
            throw DREAMException(
                "Knock-on orbit quadrature produced numerator %.17g at ir=%ld, "
                "u=%.17g, chi=%.17g, cell=[%.17g,%.17g].",
                integral, (long)ir, primaryPitch, xiStar, piece.cell.lower, piece.cell.upper
            );
        }
        numerators[piece.row] += piece.weight * integral;
    }
}

// A target piece describes the outgoing geometry before primary-pitch
// weights and trapped multiplicities are applied. Its width is measured
// relative to the receiving row, including for reassigned boundary pieces.
struct TargetPiece {
    Quadrature::Interval cell;
    len_t row;
    real_t rowWidthFraction;
    bool trapped;
};

void AddTargetPieces(
    std::vector<TargetPiece> &pieces, Quadrature::Interval cell, len_t row,
    real_t rowWidth, real_t xiT
) {
    const auto cuts = Quadrature::SplitPoints(cell, {-xiT, 0, xiT});
    for (len_t k = 1; k < cuts.size(); ++k) {
        const real_t lower = cuts[k - 1];
        const real_t upper = cuts[k];
        if (lower >= -xiT && upper <= 0)
            continue;

        const bool trapped = lower >= 0 && upper <= xiT;
        pieces.push_back({{lower, upper}, row, (upper - lower) / rowWidth, trapped});
    }
}

/**
 * Integrate A_jl over the represented part of the primary cell.
 * Use an outer partition for each outgoing piece, then evaluate the inner
 * poloidal integrals together for pieces that share a primary-pitch node.
 */
void IntegratePrimaryCell(
    len_t ir, real_t xiStar, Quadrature::Interval primary, const std::vector<TargetPiece> &targets,
    const FVM::RadialGrid *rg, Quadrature::OrbitGeometry &geometry,
    const Quadrature::GaussLegendreRule &poloidalRule,
    const Quadrature::GaussLegendreRule &pitchRule, real_t *numerators
) {
    const real_t betaMax = geometry.GetBetaMax();
    const real_t xiT = std::sqrt(betaMax);

    std::map<real_t, std::vector<OrbitPiece>> piecesByPitch;
    for (const auto &target : targets) {
        const auto cuts = Quadrature::PrimaryPitchSplitPoints(
            primary, xiStar, betaMax, {target.cell.lower, target.cell.upper}
        );
        for (len_t k = 1; k < cuts.size(); ++k) {
            const Quadrature::Interval panel{cuts[k - 1], cuts[k]};
            // Positive trapped pitch represents both orbit directions.
            if (panel.lower >= -xiT && panel.upper <= 0)
                continue;

            for (const auto &point : pitchRule.MappedPoints(panel)) {
                const real_t primaryPitch = point.x;
                // A panel near the separatrix can round onto its singular
                // pointwise orbit measure. Skip such unresolved nodes.
                if (betaMax > 0 && (std::abs(primaryPitch) == xiT
                    || primaryPitch * primaryPitch == betaMax)) {
                    continue;
                }

                const bool primaryTrapped = primaryPitch >= 0 && primaryPitch < xiT;
                const int multiplicity = primaryTrapped && target.trapped ? 2 : 1;
                piecesByPitch[primaryPitch].push_back({
                    target.cell, target.row, point.weight * target.rowWidthFraction * multiplicity,
                    primaryTrapped || target.trapped
                });
            }
        }
    }

    for (const auto &sample : piecesByPitch) {
        IntegratePieces(
            ir, xiStar, sample.first, sample.second, rg, geometry, poloidalRule, numerators
        );
    }
}

/**
 * Assemble outgoing rows of the raw cell-averaged kernel A_jl / Hgrid_l.
 * The discrete primary cell measure Hgrid_l is applied after integration,
 * so the integrand can carry the full orbit measure near the separatrix.
 */
void EvaluateGridRows(
    len_t ir, real_t xiStar, len_t l, const FVM::Grid *gridK, const FVM::Grid *gridP,
    len_t firstRow, len_t endRow, Quadrature::OrbitGeometry &geometry,
    const Quadrature::GaussLegendreRule &poloidalRule, const Quadrature::GaussLegendreRule &pitchRule,
    real_t *out
) {
    const len_t nRows = endRow - firstRow;
    std::fill(out, out + nRows, 0);
    const auto *mgP = gridP->GetMomentumGrid(ir);
    const auto *mgK = gridK->GetMomentumGrid(ir);
    const real_t primaryMetric = gridP->GetVpOverP2AtZero(ir)[l];

    // A void primary column carries no stored particles, even where its
    // interval includes passing pitch.
    if (primaryMetric == 0)
        return;

    const real_t Hgrid = mgP->GetDp2(l) * primaryMetric;
    if (!(Hgrid > 0) || !std::isfinite(Hgrid))
        throw DREAMException("Knock-on primary cell has invalid volume.");

    const real_t xiT = std::sqrt(geometry.GetBetaMax());
    std::vector<TargetPiece> targets;
    for (len_t j = firstRow; j < endRow; ++j) {
        if (gridK->GetVpOverP2AtZero(ir)[j] == 0)
            continue;

        const real_t lower = mgK->GetP2_f(j);
        const real_t upper = mgK->GetP2_f(j + 1);
        const real_t rowWidth = mgK->GetDp2(j);
        AddTargetPieces(targets, {lower, upper}, j - firstRow, rowWidth, xiT);

        // Fold the passing part of the void row crossing -xiT into the row
        // crossing +xiT. The receiving width preserves the integrated
        // outgoing source when those two rows have different widths.
        if (xiT > 0 && lower <= xiT && xiT < upper) {
            for (len_t k = 0; k < mgK->GetNp2(); ++k) {
                if (mgK->GetP2_f(k) < -xiT && mgK->GetP2_f(k + 1) > -xiT
                    && gridK->GetVpOverP2AtZero(ir)[k] == 0) {
                    AddTargetPieces(
                        targets, {mgK->GetP2_f(k), -xiT}, j - firstRow, rowWidth, xiT
                    );
                }
            }
        }
    }

    const Quadrature::Interval primary{mgP->GetP2_f(l), mgP->GetP2_f(l + 1)};
    IntegratePrimaryCell(
        ir, xiStar, primary, targets, gridP->GetRadialGrid(), geometry,
        poloidalRule, pitchRule, out
    );
    for (len_t j = 0; j < nRows; ++j) {
        out[j] /= Hgrid;
    }
}

}  // namespace

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
    Quadrature::OrbitGeometry geometry(grid_primary->GetRadialGrid()->GetFluxSurfaceAverager(), ir);
    const Quadrature::GaussLegendreRule poloidalRule(poloidalOrder);
    const Quadrature::GaussLegendreRule pitchRule(pitchOrder);
    real_t result;
    EvaluateGridRows(
        ir, xi_star, l, grid_knockon, grid_primary, j, j + 1, geometry, poloidalRule, pitchRule,
        &result
    );
    return result;
}

real_t KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, len_t poloidalOrder, len_t pitchOrder
) {
    Quadrature::OrbitGeometry geometry(grid_primary->GetRadialGrid()->GetFluxSurfaceAverager(), ir);
    const Quadrature::GaussLegendreRule poloidalRule(poloidalOrder);
    const Quadrature::GaussLegendreRule pitchRule(pitchOrder);
    return SetDeltaMatrixColumnOnGrid(
        ir, xi_star, l, grid_knockon, grid_primary, deltaCol, geometry, poloidalRule, pitchRule
    );
}

/**
 * Normalize a complete angular column to conserve the outgoing particle
 * count. Return its width-weighted integral Hquad_l / Hgrid_l so the raw
 * quadrature measure can be checked independently.
 */
real_t KnockOnUtilities::SetDeltaMatrixColumnOnGrid(
    len_t ir, real_t xi_star, len_t l, const FVM::Grid *grid_knockon, const FVM::Grid *grid_primary,
    real_t *deltaCol, Quadrature::OrbitGeometry &geometry,
    const Quadrature::GaussLegendreRule &poloidalRule, const Quadrature::GaussLegendreRule &pitchRule
) {
    const auto *mg = grid_knockon->GetMomentumGrid(ir);
    const len_t nRows = mg->GetNp2();
    EvaluateGridRows(
        ir, xi_star, l, grid_knockon, grid_primary, 0, nRows, geometry, poloidalRule, pitchRule,
        deltaCol
    );
    if (grid_primary->GetVpOverP2AtZero(ir)[l] == 0)
        return 0;

    real_t columnIntegral = 0;
    for (len_t j = 0; j < nRows; ++j) {
        columnIntegral += mg->GetDp2(j) * deltaCol[j];
    }
    if (!(columnIntegral > 0) || !std::isfinite(columnIntegral)) {
        throw DREAMException(
            "Knock-on quadrature expected a positive column at ir=%ld, l=%ld, xi_star=%.16g.",
            (long)ir, (long)l, xi_star
        );
    }

    for (len_t j = 0; j < nRows; ++j) {
        deltaCol[j] /= columnIntegral;
    }
    return columnIntegral;
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
