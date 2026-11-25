/** Regression tests for the finite-detrapping RR closure and its nonlinear coupling. */
#include <algorithm>
#include <cmath>
#include <memory>
#include <limits>
#include <vector>
#include "TrappingLimitedRRTransport.hpp"
#include "DREAM/Equations/Kinetic/TrappingLimitedRRTransport.hpp"
#include "DREAM/Equations/TransportBC.hpp"
#include "FVM/Matrix.hpp"

using namespace DREAM;
using namespace DREAMTESTS::_DREAM;

bool DREAMTESTS::_DREAM::TrappingLimitedRRTransport::Run(bool) {
    // Independent tabulated trapped-cone Sturm--Liouville eigenvalues.
    const real_t reference[][2] = {
        {0.050,492.8187397712}, {0.100,122.7078307006}, {0.221,24.5933569533},
        {0.300,13.0372539118}, {0.374,8.1438755070}, {0.472,4.8515970130},
        {0.547,3.4269917244}, {0.609,2.6196414246}, {0.662,2.0973190285},
        {0.750,1.4514192343}, {0.850,0.9232717823}
    };
    for (const auto& pair : reference) {
        const real_t xi = pair[0], expected = pair[1];
        const real_t value = DREAM::TrappingLimitedRRTransport::DetrappingEigenvalue(xi);
        const real_t hi = M_PI*M_PI/(8*xi*xi), lo = (1-xi*xi)*hi;
        if (std::abs(value-expected) > 1e-9*expected || value < lo || value > hi) {
            PrintError("Incorrect trapped-cone eigenvalue: xiT=%.4f kappa=%.13g expected=%.13g.",xi,value,expected);
            return false;
        }
    }
    for (real_t invalid : {0.98, 1.0, -0.1, std::numeric_limits<real_t>::infinity(),
                           std::numeric_limits<real_t>::quiet_NaN()}) {
        bool rejected = false;
        try { DREAM::TrappingLimitedRRTransport::DetrappingEigenvalue(invalid); }
        catch (const DREAMException&) { rejected = true; }
        if (!rejected) { PrintError("Invalid trapped boundary was accepted."); return false; }
    }
    if (DREAM::TrappingLimitedRRTransport::DetrappingEigenvalue(0) != 0)
        return false;
    std::unique_ptr<DREAM::FVM::Grid> cylindrical(InitializeGridRCylPXi(3, 6, 1, 2, 0.1, 2));
    std::unique_ptr<DREAM::FVM::Grid> toroidal(InitializeGridGeneralRPXi(3, 6, 1, 80, 20, 0.1, 2));
    bool success = CheckGrid(cylindrical.get());
    success = CheckGrid(toroidal.get()) && success;
    if (success)
        PrintOK("RR limits, finite-difference Jacobians, conservation and grid rebuild passed.");
    return success;
}

bool DREAMTESTS::_DREAM::TrappingLimitedRRTransport::CheckGrid(DREAM::FVM::Grid *grid) {
    const len_t nr = grid->GetNr(), np = grid->GetMomentumGrid(0)->GetNp1();
    std::unique_ptr<DREAM::FVM::Grid> fluid(InitializeFluidGrid(nr));
    DREAM::FVM::UnknownQuantityHandler unknowns;
    const len_t id_n = unknowns.InsertUnknown(OptionConstants::UQTY_N_COLD, "", fluid.get());
    const len_t id_T = unknowns.InsertUnknown(OptionConstants::UQTY_T_COLD, "", fluid.get());
    const len_t id_i = unknowns.InsertUnknown(OptionConstants::UQTY_ION_SPECIES, "", fluid.get(), 2);
    unknowns.InsertUnknown(OptionConstants::UQTY_F_HOT, "", grid);
    std::vector<real_t> ncold(nr), temperature(nr), ions(2*nr);
    for (len_t ir = 0; ir < nr; ir++) {
        ncold[ir] = (1+ir)*1e19;
        temperature[ir] = 100*(1+ir);
        ions[ir] = 1e18;
        ions[nr+ir] = ncold[ir];
    }
    unknowns.SetInitialValue(id_n, ncold.data());
    unknowns.SetInitialValue(id_T, temperature.data());
    unknowns.SetInitialValue(id_i, ions.data());
    const len_t *Z = new len_t[1]{1};
    std::vector<std::string> names = {"D"}, empty;
    IonHandler ionHandler(fluid->GetRadialGrid(), &unknowns, Z, 1, names, empty, empty);
    ionHandler.Rebuild();
    CollisionQuantityHandler collisions(grid, &unknowns, &ionHandler,
        OptionConstants::MOMENTUMGRID_TYPE_PXI, new CollisionQuantity::collqty_settings);
    collisions.Rebuild();

    auto perturbation = [nr](real_t amplitude) {
        auto *data = new real_t[nr+1];
        std::fill(data, data+nr+1, amplitude);
        return new DREAM::FVM::Interpolator1D(1, nr+1, new real_t[1]{0}, data);
    };
    RechesterRosenbluthTransport rr(grid, OptionConstants::MOMENTUMGRID_TYPE_PXI, perturbation(0.01));
    rr.AllocateCoefficients();
    rr.Rebuild(0, 1, &unknowns);

    DREAM::TrappingLimitedRRTransport term(grid, OptionConstants::MOMENTUMGRID_TYPE_PXI,
        perturbation(0.01), &collisions, &unknowns);
    term.AllocateCoefficients();
    TransportDiffusiveBC boundary(grid, &term, TransportDiffusiveBC::TRANSPORT_BC_F0);
    auto rebuild = [&]() {
        ionHandler.Rebuild();
        collisions.Rebuild();
        term.ResetCoefficients();
        term.Rebuild(0, 1, &unknowns);
    };
    rebuild();
    const auto *rg = grid->GetRadialGrid();
    const real_t *w = term.GetRadialJacobianInterpolationCoeffs();
    real_t *const* nu = collisions.GetNuD()->GetValue_f2();
    const real_t a = rg->GetMinorRadius();
    for (len_t ir = 0; ir <= nr; ir++) {
        const real_t xi = rg->GetXi0TrappedBoundary_fr(ir);
        const real_t kappa = DREAM::TrappingLimitedRRTransport::DetrappingEigenvalue(xi);
        if (term.GetDetrappingEigenvalue()[ir] != kappa) {
            PrintError("Incorrect trapped-cone eigenvalue at face %zu.", ir);
            return false;
        }
        for (len_t i = 0; i < np; i++) {
            const real_t faceNu = (ir < nr ? w[ir]*nu[ir][i] : 0) +
                (ir > 0 ? (1-w[ir])*nu[ir-1][i] : 0);
            const real_t cap = a*a*kappa*faceNu/5.783185962946785;
            const real_t Drr = rr.GetDiffusionCoeffRR(ir)[i];
            const real_t expected = xi == 0 ? Drr : Drr/(1+Drr/cap);
            const real_t D = term.GetDiffusionCoeffRR(ir)[i];
            if (!std::isfinite(D) || std::abs(D-expected) > 1e-12*Drr ||
                std::abs(term.GetDetrappingLimit()[ir][i]-cap) > 1e-12*cap) {
                PrintError("Incorrect eigenmode diffusion at face %zu, momentum %zu.", ir, i);
                return false;
            }
        }
    }

    std::vector<real_t> f(nr*np), residual(nr*np, 0), magnitude(nr*np);
    for (len_t ir = 0; ir < nr; ir++)
        for (len_t i = 0; i < np; i++)
            f[ir*np+i] = (1+ir*ir)*std::exp(-0.2*i);
    term.SetVectorElements(residual.data(), f.data());
    std::transform(residual.begin(), residual.end(), magnitude.begin(), [](real_t v) { return std::abs(v); });
    if (std::abs(grid->Integral(residual.data())) > 1e-12*grid->Integral(magnitude.data())) {
        PrintError("Closed-boundary radial diffusion does not conserve particles.");
        return false;
    }

    // Perturb each fluid unknown independently, including radial gradients
    // and both charge states, and compare the complete residual Jacobian.
    for (len_t id : {id_n, id_T, id_i}) {
        len_t columns = unknowns.GetUnknown(id)->NumberOfElements();
        // FVM's matrix constructor initializes all diagonal entries,
        // so use a square container for this rectangular block.
        DREAM::FVM::Matrix jac(nr*np, nr*np, nr*np);
        term.SetJacobianBlock(unknowns.GetUnknownID(OptionConstants::UQTY_F_HOT), id, &jac, f.data());
        boundary.AddToJacobianBlock(unknowns.GetUnknownID(OptionConstants::UQTY_F_HOT), id, &jac, f.data());
        jac.Assemble();
        std::vector<real_t> original(unknowns.GetUnknownData(id), unknowns.GetUnknownData(id)+columns);
        for (len_t col = 0; col < columns; col++) {
            real_t step = original[col]*1e-4;
            auto changed = original;
            changed[col] += step;
            unknowns.SetInitialValue(id, changed.data());
            rebuild();
            std::vector<real_t> plus(nr*np, 0), minus(nr*np, 0);
            term.SetVectorElements(plus.data(), f.data());
            boundary.AddToVectorElements(plus.data(), f.data());
            changed[col] -= 2*step;
            unknowns.SetInitialValue(id, changed.data());
            rebuild();
            term.SetVectorElements(minus.data(), f.data());
            boundary.AddToVectorElements(minus.data(), f.data());
            real_t error = 0, scale = 0;
            for (len_t row = 0; row < nr*np; row++) {
                real_t numeric = (plus[row]-minus[row])/(2*step);
                real_t analytic = jac.GetElement(row, col);
                error = std::max(error, std::abs(numeric-analytic));
                scale = std::max(scale, std::max(std::abs(numeric), std::abs(analytic)));
            }
            if (error > 2e-5*scale + 1e-30) {
                PrintError("Jacobian mismatch: unknown=%zu column=%zu relative error=%.3e.", id, col, error/scale);
                return false;
            }
        }
        unknowns.SetInitialValue(id, original.data());
        rebuild();
    }
    const real_t before = term.GetDiffusionCoeffRR(1)[1];
    term.GridRebuilt();
    rebuild();
    if (std::abs(before-term.GetDiffusionCoeffRR(1)[1]) > 1e-12*before) {
        PrintError("Grid rebuild changed the coefficient.");
        return false;
    }

    DREAM::TrappingLimitedRRTransport zero(grid, OptionConstants::MOMENTUMGRID_TYPE_PXI,
        perturbation(0), &collisions, &unknowns);
    zero.AllocateCoefficients();
    zero.Rebuild(0, 1, &unknowns);
    for (len_t ir = 0; ir <= nr; ir++)
        for (len_t i = 0; i < np; i++)
            if (zero.GetDiffusionCoeffRR(ir)[i] != 0 ||
                zero.GetDetrappingLimit()[ir][i] != term.GetDetrappingLimit()[ir][i]) {
                PrintError("Zero magnetic perturbation must preserve the collisional cap.");
                return false;
            }

    // Isolate the scattering-frequency limit at the outer face.
    if (rg->GetXi0TrappedBoundary_fr(nr) > 0) {
        const len_t i = np/2, last = nr-1, row = last*np+i;
        const real_t original = nu[last][i];
        const real_t capPerNu = term.GetDetrappingLimit()[nr][i] / original;
        const real_t step = 1e-7 * rr.GetDiffusionCoeffRR(nr)[i] / capPerNu;
        auto setFrequency = [&](real_t frequency) {
            nu[last][i] = frequency;
            term.ResetCoefficients();
            term.Rebuild(0, 1, &unknowns);
        };
        for (real_t frequency : {-step, 0.0}) {
            setFrequency(frequency);
            if (term.GetDiffusionCoeffRR(nr)[i] != 0 || term.GetDetrappingLimit()[nr][i] != 0) {
                PrintError("Nonpositive scattering gave a finite trapped-face flux.");
                return false;
            }
        }
        DREAM::FVM::Matrix jac(nr*np, nr*np, nr*np);
        boundary.AddToJacobianBlock(unknowns.GetUnknownID(OptionConstants::UQTY_F_HOT), id_n, &jac, f.data());
        jac.Assemble();
        const real_t analytic = jac.GetElement(row, last);
        const len_t stride = grid->GetMomentumGrid(last)->GetNCells_f2();
        const real_t dNuDn = collisions.GetNuD()->GetUnknownPartialContribution(id_n, FVM::FLUXGRIDTYPE_P2)[last*stride+i];
        setFrequency(step);
        std::vector<real_t> flux(nr*np, 0);
        boundary.AddToVectorElements(flux.data(), f.data());
        const real_t numeric = flux[row] * dNuDn / step;
        if (numeric == 0 || std::abs(analytic-numeric) > 1e-6*std::abs(numeric)) {
            PrintError("Incorrect positive-side boundary Jacobian at zero scattering.");
            return false;
        }
        setFrequency(original);
    }
    return true;
}
