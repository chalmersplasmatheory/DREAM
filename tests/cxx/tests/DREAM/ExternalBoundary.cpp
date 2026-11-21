#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>
#include "ExternalBoundary.hpp"
#include "DREAM/Equations/Kinetic/ExternalKineticBoundary.hpp"
#include "FVM/BlockMatrix.hpp"
#include "FVM/Grid/AnalyticBRadialGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
#include "FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
#include "FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
#include "FVM/Equation/AdvectionTerm.hpp"
#include "FVM/Equation/DiffusionTerm.hpp"

using namespace DREAM;
using namespace DREAMTESTS::_DREAM;

namespace {
using Boundary = ExternalKineticBoundary;

struct GridDeleter {
    void operator()(FVM::Grid *grid) const {
        auto *radial = grid->GetRadialGrid();
        delete grid;
        delete radial;
    }
};
using OwnedGrid = std::unique_ptr<FVM::Grid, GridDeleter>;

OwnedGrid ToroidalGrid() {
    using namespace FVM;
    static const real_t r[] = {0, 1}, one[] = {1, 1}, zero[] = {0, 0};
    auto *profiles = new AnalyticBRadialGridGenerator::shape_profiles{
        2, 2, 2, 2, 2, one, r, zero, r, one, r, zero, r, zero, r
    };
    auto *radial = new RadialGrid(new AnalyticBRadialGridGenerator(2, 0.19, 0.21, 1, 128, profiles));
    auto *momentum = new PXiGrid::PXiMomentumGrid(new PXiGrid::MomentumGridGenerator(
        new PXiGrid::PUniformGridGenerator(4, 0, 2), new PXiGrid::XiUniformGridGenerator(9)
    ), 0, radial);
    OwnedGrid grid(new Grid(radial, momentum));
    grid->Rebuild(0);
    return grid;
}

// Depend on a fluid unknown to exercise the coefficient Jacobian, including
// the derivative of the limiter's total inward-rate denominator.
class BoundaryAdvection : public FVM::AdvectionTerm {
    len_t id;
public:
    BoundaryAdvection(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t id)
        : AdvectionTerm(g), id(id) {
        AddUnknownForJacobian(u, id);
        AllocateDifferentiationCoefficients();
    }
    void Rebuild(real_t, real_t, FVM::UnknownQuantityHandler *u) override {
        for (len_t ir = 0; ir < nr; ++ir)
            for (len_t j = 0; j < n2[ir]; ++j)
                F1(ir, n1[ir], j) = u->GetUnknownData(id)[ir] * (j%2 ? -2.0 : 1.0);
    }
    void SetPartialAdvectionTerm(len_t, len_t) override {
        ResetDifferentiationCoefficients();
        for (len_t ir = 0; ir < nr; ++ir)
            for (len_t j = 0; j < n2[ir]; ++j)
                dF1(ir, n1[ir], j, 0) = j%2 ? -2.0 : 1.0;
    }
};

class BoundaryDiffusion : public FVM::DiffusionTerm {
    len_t id;
public:
    BoundaryDiffusion(FVM::Grid *g, FVM::UnknownQuantityHandler *u, len_t id)
        : DiffusionTerm(g), id(id) {
        AddUnknownForJacobian(u, id);
        AllocateDifferentiationCoefficients();
    }
    void Rebuild(real_t, real_t, FVM::UnknownQuantityHandler *u) override {
        for (len_t ir = 0; ir < nr; ++ir)
            for (len_t j = 0; j < n2[ir]; ++j) {
                D11(ir, n1[ir], j) = 0.1*u->GetUnknownData(id)[ir];
                D12(ir, n1[ir], j) = 0.02*u->GetUnknownData(id)[ir];
            }
    }
    void SetPartialDiffusionTerm(len_t, len_t) override {
        ResetDifferentiationCoefficients();
        for (len_t ir = 0; ir < nr; ++ir)
            for (len_t j = 0; j < n2[ir]; ++j) {
                dD11(ir, n1[ir], j) = 0.1;
                dD12(ir, n1[ir], j) = 0.02;
            }
    }
};
}

bool ExternalBoundary::CheckGrid(FVM::Grid *kineticGrid, FVM::Grid *fluidGrid) {
    const len_t nr = kineticGrid->GetNr(), cells = kineticGrid->GetNCells();
    FVM::UnknownQuantityHandler unknowns;
    const len_t idF = unknowns.InsertUnknown("f_test", "Test distribution", kineticGrid);
    const len_t idP = unknowns.InsertUnknown("n_external_pos", "Positive inventory", fluidGrid);
    const len_t idN = unknowns.InsertUnknown("n_external_neg", "Negative inventory", fluidGrid);
    const len_t idT = unknowns.InsertUnknown("n_external", "Total inventory", fluidGrid);
    const len_t idA = unknowns.InsertUnknown("a_test", "Advection amplitude", fluidGrid);
    const len_t idD = unknowns.InsertUnknown("d_test", "Diffusion amplitude", fluidGrid);
    std::vector<real_t> f(cells), positive(nr, 5), negative(nr, 2), a(nr, 1), d(nr, 1);
    for (len_t ir = 0, offset = 0; ir < nr; ++ir) {
        const auto *mg = kineticGrid->GetMomentumGrid(ir);
        for (len_t j = 0; j < mg->GetNp2(); ++j)
            for (len_t i = 0; i < mg->GetNp1(); ++i)
                // Exercise all three reconstruction branches: decreasing,
                // increasing (slope set to zero), and trace clipped to zero.
                f[offset+j*mg->GetNp1()+i] = (1 + 0.3*ir + 0.1*j)
                    * std::exp((j%3 == 0 ? -0.6 : (j%3 == 1 ? 0.2 : -6.0))*mg->GetP1(i));
        offset += mg->GetNCells();
    }
    unknowns.SetInitialValue(idF, f.data());
    unknowns.SetInitialValue(idP, positive.data());
    unknowns.SetInitialValue(idN, negative.data());
    unknowns.SetInitialValue(idT, nullptr);
    unknowns.SetInitialValue(idA, a.data());
    unknowns.SetInitialValue(idD, d.data());
    unknowns.SaveStep(0, false);
    FVM::Operator op(kineticGrid);
    op.AddTerm(new BoundaryAdvection(kineticGrid, &unknowns, idA));
    op.AddTerm(new BoundaryDiffusion(kineticGrid, &unknowns, idD));
    const std::array poles{idP, idN};
    Boundary kinetic(kineticGrid, kineticGrid, &op, idF, poles, &unknowns);
    Boundary fluid(fluidGrid, kineticGrid, &op, idF, poles, &unknowns, Boundary::DENSITY);
    Boundary negativeSide(fluidGrid, kineticGrid, &op, idF, poles, &unknowns, Boundary::NEGATIVE_DENSITY);
    Boundary positiveSide(fluidGrid, kineticGrid, &op, idF, poles, &unknowns, Boundary::POSITIVE_DENSITY);
    const real_t dt = 0.3;
    auto rebuild = [&]() {
        op.RebuildTerms(dt, dt, &unknowns);
        kinetic.Rebuild(dt, &unknowns);
        fluid.Rebuild(dt, &unknowns);
        negativeSide.Rebuild(dt, &unknowns);
        positiveSide.Rebuild(dt, &unknowns);
    };
    rebuild();
    std::vector<real_t> rate(nr);
    kinetic.GetRequestedInflow(rate.data());
    std::vector<real_t> k(cells), balance(nr), positiveBalance(nr), negativeBalance(nr), integral(nr), net(nr), out(nr), in(nr);
    std::array inward{std::vector<real_t>(nr), std::vector<real_t>(nr)};
    auto requested = inward, netByPole = inward;
    for (auto pole : {Boundary::POSITIVE, Boundary::NEGATIVE}) {
        kinetic.GetRequestedInflow(requested[pole].data(), pole);
        kinetic.GetFlux(netByPole[pole].data(), nullptr, inward[pole].data(), pole);
    }
    kinetic.AddToVectorElements(k.data(), f.data());
    fluid.AddToVectorElements(balance.data(), f.data());
    negativeSide.AddToVectorElements(negativeBalance.data(), f.data());
    positiveSide.AddToVectorElements(positiveBalance.data(), f.data());
    kineticGrid->IntegralMomentum(k.data(), integral.data());
    kinetic.GetFlux(net.data(), out.data(), in.data());
    for (len_t ir = 0; ir < nr; ++ir) {
        const real_t scale = std::max(real_t(1), out[ir]+in[ir]);
        if (std::abs(integral[ir]+balance[ir]) > 1e-12*scale
            || std::abs(net[ir]-integral[ir]) > 1e-12*scale
            || std::abs(net[ir]-netByPole[0][ir]-netByPole[1][ir]) > 1e-12*scale
            || std::abs(positiveBalance[ir]+netByPole[0][ir]) > 1e-12*scale
            || std::abs(negativeBalance[ir]+netByPole[1][ir]) > 1e-12*scale
            || std::abs(in[ir]-inward[0][ir]-inward[1][ir]) > 1e-12*scale) {
            PrintError("The total and tagged external balances must use the kinetic flux.");
            return false;
        }
        const std::array density{positive[ir], negative[ir]};
        for (auto pole : {Boundary::POSITIVE, Boundary::NEGATIVE}) {
            const real_t I = requested[pole][ir], n = density[pole];
            if (I <= 0 || inward[pole][ir]*dt > n
                || std::abs(inward[pole][ir] - I*n/(n+dt*I)) > 1e-12*scale) {
                PrintError("Inflow must be limited by the inventory of its own pole.");
                return false;
            }
        }
    }

    // Allocate each block uniformly, as in the solver.
    FVM::BlockMatrix jac;
    const len_t kineticBlock = jac.CreateSubEquation(cells, 1+kinetic.GetNumberOfNonZerosPerRow_jac(), idF);
    const len_t positiveBlock = jac.CreateSubEquation(nr, 1+positiveSide.GetNumberOfNonZerosPerRow_jac(), idP);
    const len_t negativeBlock = jac.CreateSubEquation(nr, 1+negativeSide.GetNumberOfNonZerosPerRow_jac(), idN);
    const len_t totalBlock = jac.CreateSubEquation(nr, 1+fluid.GetNumberOfNonZerosPerRow_jac(), idT);
    const len_t advectionBlock = jac.CreateSubEquation(nr, 1, idA);
    const len_t diffusionBlock = jac.CreateSubEquation(nr, 1, idD);
    const std::array blocks{
        std::pair{idF, kineticBlock}, std::pair{idP, positiveBlock}, std::pair{idN, negativeBlock},
        std::pair{idT, totalBlock}, std::pair{idA, advectionBlock}, std::pair{idD, diffusionBlock}
    };
    jac.ConstructSystem();

    auto assembleBoundaryJacobian = [&]() {
        for (const auto &[id, block] : blocks) {
            jac.SelectSubEquation(kineticBlock, block);
            kinetic.AddToJacobianBlock(idF, id, &jac, f.data());
            jac.SelectSubEquation(positiveBlock, block);
            positiveSide.AddToJacobianBlock(idP, id, &jac, f.data());
            jac.SelectSubEquation(negativeBlock, block);
            negativeSide.AddToJacobianBlock(idN, id, &jac, f.data());
            jac.SelectSubEquation(totalBlock, block);
            fluid.AddToJacobianBlock(idT, id, &jac, f.data());
        }
        jac.Assemble();
    };

    // Current iterates change outflow, while inflow uses the accepted trace.
    const std::vector<real_t> zero(nr, 0);
    auto iterate = f;
    for (auto &value : iterate)
        value *= 2;
    unknowns.Store(idF, iterate.data());
    unknowns.Store(idD, zero.data());
    rebuild();
    kinetic.GetRequestedInflow(rate.data());
    for (len_t ir = 0, offset = 0; ir < nr; ++ir) {
        const auto *mg = kineticGrid->GetMomentumGrid(ir);
        const len_t np = mg->GetNp1();
        const auto *advection = op.GetAdvectionCoeff1(ir);
        real_t expected = 0;
        for (len_t j = 0; j < mg->GetNp2(); ++j) {
            const len_t cell = j*np+np-1, face = j*(np+1)+np;
            if (kineticGrid->IsNegativePitchTrappedIgnorableCell(ir, j)
                || kineticGrid->GetVp(ir)[cell] <= 0)
                continue;
            const real_t weight = kineticGrid->GetVp_f1(ir)[face]*mg->GetDp2(j)/kineticGrid->GetVpVol(ir);
            expected += weight*std::max(real_t(0), -advection[face])*f[offset+cell];
        }
        if (std::abs(rate[ir]-expected) > 1e-12*std::max(real_t(1), expected)) {
            PrintError("External boundary inward advection does not use the accepted last cell.");
            return false;
        }
        offset += mg->GetNCells();
    }
    const auto acceptedRate = rate;
    kinetic.GetFlux(net.data(), out.data(), in.data());
    for (len_t ir=0, offset=0; ir<nr; ++ir) {
        const auto *mg = kineticGrid->GetMomentumGrid(ir);
        const len_t np = mg->GetNp1();
        const auto *advection = op.GetAdvectionCoeff1(ir);
        real_t expected = 0;
        for (len_t j=0; j<mg->GetNp2(); ++j) {
            if (kineticGrid->IsNegativePitchTrappedIgnorableCell(ir, j)
                || kineticGrid->GetVp(ir)[j*np+np-1] <= 0)
                continue;
            const len_t cell = offset+j*np+np-1, face = j*(np+1)+np;
            const real_t weight = kineticGrid->GetVp_f1(ir)[face]*mg->GetDp2(j)/kineticGrid->GetVpVol(ir);
            expected += weight*std::max(real_t(0), advection[face])*iterate[cell];
        }
        if (std::abs(out[ir]-expected) > 1e-12*std::max(real_t(1), expected)) {
            PrintError("External boundary outflow lost its current-distribution dependence.");
            return false;
        }
        offset += mg->GetNCells();
    }

    auto negativeIterate = iterate;
    for (auto &value : negativeIterate)
        value = -value;
    unknowns.Store(idF, negativeIterate.data());
    rebuild();
    kinetic.GetRequestedInflow(rate.data());
    kinetic.GetFlux(net.data(), out.data(), in.data());
    for (len_t ir=0; ir<nr; ++ir)
        if (std::abs(rate[ir]-acceptedRate[ir]) > 1e-12*std::max(real_t(1), acceptedRate[ir])
            || out[ir] != 0 || in[ir] <= 0) {
            PrintError("A negative current iterate altered the accepted inflow request.");
            return false;
        }

    // Advancing and rolling back a step must advance and restore the trace.
    unknowns.Store(idF, iterate.data());
    unknowns.SaveStep(dt, false);
    // Flux output is written after SaveStep, before the next Rebuild.
    kinetic.GetRequestedInflow(rate.data());
    for (len_t ir=0; ir<nr; ++ir)
        if (std::abs(rate[ir]-acceptedRate[ir]) > 1e-12*std::max(real_t(1), acceptedRate[ir]))
            return false;
    op.RebuildTerms(2*dt, dt, &unknowns);
    kinetic.Rebuild(2*dt, &unknowns);
    kinetic.GetRequestedInflow(rate.data());
    for (len_t ir=0; ir<nr; ++ir)
        if (std::abs(rate[ir]-2*acceptedRate[ir]) > 1e-12*std::max(real_t(1), rate[ir]))
            return false;
    unknowns.RollbackSaveStep();
    rebuild();
    kinetic.GetRequestedInflow(rate.data());
    for (len_t ir=0; ir<nr; ++ir)
        if (std::abs(rate[ir]-acceptedRate[ir]) > 1e-12*std::max(real_t(1), acceptedRate[ir]))
            return false;
    // Exercise all Jacobian blocks with different current and accepted traces.
    f = iterate;
    unknowns.Store(idF, f.data());

    // Start from zero face fluxes, then use the same matrix with both signs.
    // The initial sparsity pattern must include derivatives activated later.
    unknowns.Store(idA, zero.data());
    rebuild();
    assembleBoundaryJacobian();
    jac.Zero();
    unknowns.Store(idA, a.data());
    unknowns.Store(idD, d.data());
    rebuild();
    assembleBoundaryJacobian();

    // Perturb every column, including remote pitch cells and other radii.
    auto residual = [&](std::vector<real_t>& values) {
        const real_t *distribution = unknowns.GetUnknownData(idF);
        kinetic.AddToVectorElements(values.data(), distribution);
        positiveSide.AddToVectorElements(values.data()+cells, distribution);
        negativeSide.AddToVectorElements(values.data()+cells+nr, distribution);
        fluid.AddToVectorElements(values.data()+cells+2*nr, distribution);
    };
    for (const auto &[id, block] : blocks) {
        const len_t columns = unknowns.GetUnknown(id)->NumberOfElements();
        const len_t columnOffset = jac.GetOffset(block);
        auto *x = unknowns.GetUnknownData(id);
        for (len_t col = 0; col < columns; ++col) {
            const real_t original = x[col], h = 1e-5*std::max(real_t(1), std::abs(original));
            std::vector<real_t> plus(cells+3*nr), minus(cells+3*nr);
            x[col] = original+h;
            rebuild();
            residual(plus);
            x[col] = original-h;
            rebuild();
            residual(minus);
            x[col] = original;
            rebuild();
            for (len_t row = 0; row < cells+3*nr; ++row) {
                const real_t actual = jac.GetElement(row, columnOffset+col);
                const real_t numerical = (plus[row]-minus[row])/(2*h);
                if (std::abs(actual-numerical) > 2e-7*std::max({real_t(1), std::abs(actual), std::abs(numerical)})) {
                    PrintError("External-boundary Jacobian mismatch: unknown %u, row %u, column %u: %g vs %g.",
                        unsigned(id), unsigned(row), unsigned(col), actual, numerical);
                    return false;
                }
            }
        }
    }

    // An empty (or slightly negative Newton-iterate) reservoir must never
    // supply inflow. Nor may a finite reservoir invent an empty trace.
    for (real_t inventory : {real_t(0), real_t(-1), real_t(1e-15), real_t(1e15)}) {
        std::fill(negative.begin(), negative.end(), inventory);
        unknowns.Store(idN, negative.data());
        rebuild();
        kinetic.GetFlux(nullptr, nullptr, in.data(), Boundary::NEGATIVE);
        kinetic.GetFlux(nullptr, nullptr, inward[0].data(), Boundary::POSITIVE);
        for (len_t ir = 0; ir < nr; ++ir)
            if (in[ir] < 0 || in[ir]*dt > std::max(real_t(0), inventory)*(1+1e-12)
                || inward[0][ir] <= 0)
                return false;
    }
    std::fill(f.begin(), f.end(), 0);
    unknowns.Store(idF, f.data());
    rebuild();
    kinetic.GetFlux(net.data(), out.data(), in.data());
    for (len_t ir = 0; ir < nr; ++ir)
        if (out[ir] != 0 || in[ir] <= 0)
            return false;
    unknowns.SaveStep(dt, false);
    op.RebuildTerms(2*dt, dt, &unknowns);
    kinetic.Rebuild(2*dt, &unknowns);
    kinetic.GetFlux(net.data(), out.data(), in.data());
    for (len_t ir = 0; ir < nr; ++ir)
        if (out[ir] != 0 || in[ir] != 0)
            return false;
    // Isolate one outward face at each radius. In the cylinder the central
    // cell straddles zero. In the torus this test selects a trapped cell.
    unknowns.Store(idD, zero.data());
    for (real_t positiveShare : {real_t(0), real_t(0.5), real_t(1)}) {
        std::fill(f.begin(), f.end(), 0);
        for (len_t ir = 0, offset = 0; ir < nr; ++ir) {
            const auto *mg = kineticGrid->GetMomentumGrid(ir);
            const real_t xiT = kineticGrid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
            const len_t np = mg->GetNp1();
            bool found = false;
            for (len_t j = 0; j < mg->GetNp2(); ++j) {
                const real_t lower = mg->GetP2_f(j), upper = mg->GetP2_f(j+1);
                if (kineticGrid->IsNegativePitchTrappedIgnorableCell(ir, j)
                    || kineticGrid->GetVp(ir)[j*np+np-1] <= 0)
                    continue;
                const bool shared = xiT == 0 ? lower < 0 && std::abs(lower+upper) < 1e-14
                    : lower >= -xiT && upper <= xiT;
                if ((positiveShare == 0 && upper <= -xiT)
                    || (positiveShare == 1 && lower >= xiT)
                    || (positiveShare == 0.5 && shared)) {
                    f[offset+j*np+np-1] = 1;
                    a[ir] = j%2 ? -1 : 1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                PrintError("Test grid has no active cell for the requested pole split.");
                return false;
            }
            offset += mg->GetNCells();
        }
        unknowns.Store(idA, a.data());
        unknowns.Store(idF, f.data());
        op.RebuildTerms(2*dt, dt, &unknowns);
        kinetic.Rebuild(2*dt, &unknowns);
        kinetic.GetFlux(nullptr, out.data());
        kinetic.GetFlux(netByPole[0].data(), nullptr, nullptr, Boundary::POSITIVE);
        kinetic.GetFlux(netByPole[1].data(), nullptr, nullptr, Boundary::NEGATIVE);
        for (len_t ir = 0; ir < nr; ++ir)
            if (out[ir] <= 0
                || std::abs(netByPole[0][ir]/out[ir] - positiveShare) > 1e-12
                || std::abs(netByPole[1][ir]/out[ir] - (1-positiveShare)) > 1e-12) {
                PrintError("The outward face did not route its flux to the expected poles.");
                return false;
            }
    }
    return true;
}

bool ExternalBoundary::Run(bool) {
    OwnedGrid cylinder(InitializeGridRCylPXi(2, 4, 7, 2, 0, 2));
    OwnedGrid fluid(InitializeFluidGrid(2));
    if (!CheckGrid(cylinder.get(), fluid.get()))
        return false;
    auto torus = ToroidalGrid();
    if (!CheckGrid(torus.get(), fluid.get()))
        return false;
    PrintOK("Two-pole boundary: signed balances, independent reservoirs, lagged trace and Jacobian (cylinder and torus).");
    return true;
}
