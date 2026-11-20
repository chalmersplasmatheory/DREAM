# DREAM FVM discretisation of the diffusion operator and the pitch-angle scattering time

This note documents (1) how DREAM discretises the momentum-space diffusion
operator in its finite-volume method (FVM), and (2) why the correct pitch-angle
scattering time associated with the stored coefficient `other/runaway/Dxx` is

    tau_xi ~ 1 / (2 * Dxx)          (order-of-magnitude: 1/Dxx)

and **not** `p^2 / Dxx`. The latter formula is correct only when `D` is the
*parallel-momentum* diffusion coefficient (units: momentum^2 / s), not the
pitch-angle (xi-space) coefficient output by DREAM.

Reference sources:
- `fvm/Equation/DiffusionTerm.set.cpp` — discrete flux assembly
- `fvm/Equation/DiffusionTerm.cpp`, `include/FVM/Equation/DiffusionTerm.hpp`
- `tests/cxx/tests/FVM/DiffusionTerm.cpp` — conservativity check
- `src/Equations/Kinetic/QuasilinearDiffusionTerm.cpp` — QL coefficient fill

---

## 1. Continuous form of the momentum-space diffusion operator

On the `(p, xi)` momentum grid (p dimensionless momentum `p = |p|/(m_e c)`,
`xi = p_par/p`), DREAM's kinetic equation contains the diffusion term in
conservative (flux) form:

```
df/dt = (1/Vp) * [ d/dp  ( Vp_f1 * ( D11 df/dp  + D12 df/dxi ) )
                 + d/dxi ( Vp_f2 * ( D21 df/dp  + D22 df/dxi ) ) ]
```

with

- `Vp  = p^2` — phase-space Jacobian at cell centres,
- `Vp_f1`, `Vp_f2` — the same Jacobian evaluated on the p- and xi-face grids,
- `D11 = D_pp`, `D22 = D_xixi`, `D12 = D_pxi`, `D21 = D_xip`.

In the pure `D22` (pitch-angle) channel and because `Vp = p^2` does not depend
on `xi`, the operator reduces exactly to

```
df/dt = d/dxi ( D22 * df/dxi ) .                       (1)
```

There is **no** residual `p^2` factor in (1): the Jacobian cancels between the
numerator `Vp_f2` and the denominator `Vp` of the flux divergence. This is the
key structural fact behind the scattering-time estimate.

## 2. DREAM's finite-volume differencing

DREAM is strictly conservative: the grid carries cell volumes
`Vp[i,j] * dp1[i] * dp2[j]` and every face flux leaving one cell enters the
neighbour. The unit test `CheckConservation` verifies that the column sums of
the discretisation matrix weighted by the cell volumes vanish to machine
precision (particle number is conserved up to boundary fluxes).

Concretely, in `DiffusionTerm.set.cpp` the xi-xi channel is assembled as a
three-point stencil (loop over cells `(i,j)` = `(p_i, xi_j)`). Each xi-face
carries the Fick flux (up to the framework's global sign convention)

    Phi_xi(i, j-1/2) = - D22(i, j) * Vp_f2[i,j] * ( f(i,j) - f(i,j-1) ) / dp2_f[j-1]

and the cell update is the flux difference divided by the cell volume
`Vp * dp1 * dp2` (exactly the reference form in
`tests/cxx/tests/FVM/EquationTerm.cpp::EvaluateFVM`). The code assembles, for
the face between cells `j-1` and `j`, the coupling coefficient

    S = D22(ir, i, j) * Vp_f2[i,j] / ( dp2[j] * dp2_f[j-1] * Vp[i,j] )

appearing as a `+S` / `-S` pair on the diagonal and the neighbour column.
The same pattern repeats for the `j`/`j+1` face with `dp2_f[j]`. The result is
a standard conservative three-point diffusion stencil. Notes:

- `D22` lives on **xi-faces** (`np2+1` values, layout `(nt, nr, nxi+1, np)` in
  the output), `D11` on **p-faces** (`np1+1`); cross coefficients `D12`/`D21`
  sit on p-/xi-faces respectively and use centred differences in the
  transverse direction.
- Face values of `f` are obtained by centred differences of the two adjacent
  cells (`dp2_f` = distance between cell centres), so the scheme is
  second-order accurate for smooth `f`.
- Boundary behaviour: at `xi = +-1` (first/last faces) no diffusion flux is
  assembled — the pitch-angle flux vanishes at the domain edge, consistent
  with the conservative form; at `p = 0` and `p = pmax` likewise.
- Trapped-particle cells (`Vp = 0` region) are skipped entirely.

Because the stencil is in flux form, any monotone `f` diffuses without
creating or destroying particles; the discrete operator is exactly (1)
discretised.

## 3. From the operator to the scattering time

Equation (1) is a 1D diffusion equation in the variable `xi` with diffusion
coefficient `D22 [s^-1]` (xi is dimensionless, so Dxx carries units of 1/time).
For a constant `D22` the mean-square displacement of `xi` is

    < (xi(t) - xi(0))^2 > = 2 * D22 * t .                (2)

This follows from the Fokker-Planck form of (1) — the standard result
`d<Δξ²>/dt = 2D` — and is equally true of the discrete stencil in the limit
where `f` is resolved by several cells.

Two equivalent ways to read off the characteristic time:

1. **Randomisation time** (O(1) change of pitch angle, `|Δxi| ~ 1`):

       tau_xi = 1 / (2 * Dxx)   ~   1 / Dxx .

2. **Discrete hop-rate picture**: in one cell-to-cell hop the particle changes
   `xi` by `dp2_f`; the hop rate is `Dxx / dp2_f^2`, so the time to randomise
   over a range `Δxi` is `Δxi^2 / (2 Dxx)` — again `~ 1/Dxx` for `Δxi ~ 1`.

Numerical check with the values from the low-amplitude run
(`figures_wk40_A4_t3_s2_i0.2_T1000`): `Dxx ~ 1e3 s^-1` gives

    tau_xi = 1 / (2e3) = 5e-4 s  <<  0.01 s ,

so a single saved frame (10 ms) already contains many scattering times —
consistent with the observed complete pitch-angle scattering between
consecutive frames. (Resonant cells have larger Dxx and scatter even faster;
the `D_pxi` cross term contributes as well.)

## 4. Why `p^2/D` appears in the literature — and where it applies

The quasilinear tensor is often written in momentum space,

    D^(n,k) = A_{n,k} * Δ Δ^T ,   Δ = ( ω - k_par v_par ,  p_perp k_par c / gamma ),

i.e. the components are diffusion coefficients of the **momentum components**
(`D_parpar`, `D_parperp`, ...), with units momentum^2 / s. If a pitch-angle
diffusion coefficient `D_xixi` is defined from parallel-momentum diffusion at
fixed `p`, then `xi = p_par/p` implies `Δxi = Δp_par / p`, hence

    <Δxi²> = 2 D_parpar t / p^2   =>   tau_xi = p^2 / (2 D_parpar) .   (3)

So (3) is correct, but its `D` is `D_parpar` (momentum-space), while the
quantity stored in `other/runaway/Dxx` is the coefficient of
`d/dxi(Dxx df/dxi)` in (1), already the pitch-angle-space coefficient. The two
are related by

    D_parpar = p^2 * Dxx ,

so (3) with `D = Dxx` counts the factor `p^2` twice. Inserting `p = 16`,
`Dxx = 1e3 s^-1`:

- wrong (momentum-space formula with pitch-angle coefficient):
  `tau = p^2/Dxx = 256/1e3 = 0.26 s` — contradicts the observed scattering
  within 0.01 s;
- correct (equations (1)-(2)):
  `tau = 1/(2 Dxx) = 5e-4 s` — consistent with the simulation.

Note that at ultra-relativistic momenta `p >> 1` a "1 s" of xi-diffusion at
`Dxx ~ 1e3` corresponds to an enormous parallel-momentum diffusion
`D_parpar = p^2 Dxx ~ 2.6e5`, which is exactly why the two conventions can be
confused: both give the same physics when the right coefficient is paired with
the right formula.

## 5. Practical consequences

- To estimate how fast a run will isotropise in pitch angle, read
  `other/runaway/Dxx` (or the cached `D_xixi` in `QuasilinearDiffusionTerm`)
  and use `tau_xi ~ 1/(2 Dxx)`. There is no `p^2` factor.
- If an analytic estimate is built from the momentum-space tensor
  `A ΔΔ^T`, first convert to the xi-space coefficient
  `Dxx = D_parpar / p^2` (or use the compact coefficient
  `D_xixi = Sigma * Lambda^2` of the paper convention directly — that one is
  already the xi-space coefficient, units s^-1), then apply `1/(2 Dxx)`.
- For time resolution of the *transient* scattering process in outputs, the
  save interval must be `<= tau_xi`; at `Dxx ~ 1e3 s^-1` even the 1 ms
  internal step is marginal and the 10 ms save step only resolves the
  pre/post-scattering states.
