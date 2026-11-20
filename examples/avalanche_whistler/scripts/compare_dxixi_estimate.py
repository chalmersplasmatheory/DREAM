#!/usr/bin/env python3
#
# Compare the analytic D_xixi estimate (compact quasilinear coefficient,
# n=1 resonance, delta(omega) ~ 1/(|v_g - v_par cos th| dk)) against the
# value computed by DREAM's QuasilinearDiffusionTerm and stored in the
# output ('other/runaway/Dxx').
#
# Estimate (single effective mode, table values by default):
#
#   D_xixi ~ 4*pi*e^2/(m_e^2 c^2 eps0 D_diel) * W * |psi_norm|^2
#            * 1/(|v_g - v_par cos th| * dk) * Lambda^2
#   Lambda = -xi + k_par v / omega
#
# Usage:
#   python compare_dxixi_estimate.py --h5 smoke_test/dxixi_check.h5
#

import argparse
import h5py
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')

# Physical constants
E_CHARGE = 1.602176634e-19   # C
M_E = 9.10938356e-31         # kg
C_LIGHT = 299792458.0        # m/s
EPS0 = 8.854187817e-12       # F/m
M_E_C2_J = M_E * C_LIGHT**2  # J


def whistler_dispersion_simple(k, theta):
    """omega = k*|k_par|*w_factor is NOT used here; print reference only."""
    raise NotImplementedError


def estimate_dxixi(W, D_diel, psi2, Lambda, vg, vpar_costheta, dk):
    """
    Analytic estimate of the compact D_xixi coefficient:

        D ~ 4 pi e^2 / (m^2 c^2 eps0 D) * W * |psi|^2
            * 1/(|vg - v_par cos th| * dk) * Lambda^2
    """
    denom = np.abs(vg - vpar_costheta)
    pref = 4.0 * np.pi * E_CHARGE**2 / (M_E**2 * C_LIGHT**2 * EPS0 * D_diel)
    return pref * W * psi2 * (1.0 / (denom * dk)) * Lambda**2


def lambda_kinematics(p, xi, k, theta_k, omega):
    """Lambda = -xi + k_par v / omega at (p, xi), k_par = k cos(theta_k)."""
    gamma = np.sqrt(1.0 + p**2)
    v = (p / gamma) * C_LIGHT
    k_par = k * np.cos(theta_k)
    return -xi + k_par * v / omega


def resonance_check(p, xi, k, theta_k, omega, n, B0):
    """Delta_omega = omega - k_par v_par - n Omega_ce / gamma."""
    gamma = np.sqrt(1.0 + p**2)
    v_par = (p / gamma) * xi * C_LIGHT
    k_par = k * np.cos(theta_k)
    omega_ce = E_CHARGE * B0 / M_E
    return omega - k_par * v_par - n * omega_ce / gamma


def wave_space_factor(k_min, k_max, th_min, th_max):
    return (2.0 * np.pi * (k_max**3 - k_min**3) / 3.0
            * (np.cos(th_min) - np.cos(th_max)))


def load_program_dxixi(h5_path, p_target, xi_target, grid='runaway'):
    """
    Read 'other/<grid>/Dxx' from a DREAM output and return the value at the
    momentum cell center nearest to (p_target, xi_target), plus metadata.
    """
    with h5py.File(h5_path, 'r') as f:
        if f'other/{grid}/Dxx' not in f:
            print(f"Error: 'other/{grid}/Dxx' not found in {h5_path}")
            return None
        Dxx = f[f'other/{grid}/Dxx'][:]        # (nt, nr, nxi+1, np)
        p = f[f'grid/{grid}/p1'][:]
        xi = f[f'grid/{grid}/p2'][:]
        t = f['grid/t'][:]

    # average xi-faces -> cell centers (nt, nxi, np)
    Dxx_c = 0.5 * (Dxx[:, 0, :-1, :] + Dxx[:, 0, 1:, :])
    nt_q = Dxx_c.shape[0]
    if len(t) == nt_q + 1:
        t = t[1:]

    ip = int(np.argmin(np.abs(p - p_target)))
    jx = int(np.argmin(np.abs(xi - xi_target)))
    return {
        't': t, 'p_grid': p, 'xi_grid': xi, 'Dxx': Dxx_c,
        'ip': ip, 'jx': jx,
        'p_actual': p[ip], 'xi_actual': xi[jx],
        'value': float(Dxx_c[-1, jx, ip]),
        'max_all': float(np.max(np.abs(Dxx_c[-1]))),
    }


def main():
    ap = argparse.ArgumentParser(description='Compare analytic vs program D_xixi')
    ap.add_argument('--h5', required=True, help='DREAM output with other/runaway/Dxx')
    ap.add_argument('--grid', default='runaway', choices=['runaway', 'hottail'])
    # phase-space point
    ap.add_argument('--p', type=float, default=20.0)
    ap.add_argument('--xi', type=float, default=0.75)
    ap.add_argument('--n', type=int, default=1, help='resonance harmonic')
    # table values (Dxixi estimate table)
    ap.add_argument('--W', type=float, default=5e-3, help='wave energy density [J/m^3]')
    ap.add_argument('--D-diel', type=float, default=35.0, help='dielectric derivative D')
    ap.add_argument('--psi2', type=float, default=0.6, help='|psi_norm|^2')
    ap.add_argument('--Lambda', type=float, default=None,
                    help='Lambda = -xi + k_par v/omega (default: computed)')
    ap.add_argument('--vg', type=float, default=9.8e7, help='group velocity [m/s]')
    ap.add_argument('--vpar-costheta', type=float, default=-1.6e8,
                    help='v_par*cos(theta_k) [m/s]')
    ap.add_argument('--dk', type=float, default=11.6, help='k_max-k_min [1/m]')
    # spectrum used by the program run
    ap.add_argument('--k-min', type=float, default=50.2)
    ap.add_argument('--k-max', type=float, default=61.8)
    ap.add_argument('--ktheta-min', type=float, default=2.38)
    ap.add_argument('--ktheta-max', type=float, default=2.47)
    ap.add_argument('--B0', type=float, default=1.4)
    ap.add_argument('--omega', type=float, default=2 * np.pi * 476e6,
                    help='wave angular frequency (rad/s)')
    args = ap.parse_args()

    k_nom = 0.5 * (args.k_min + args.k_max)
    th_nom = 0.5 * (args.ktheta_min + args.ktheta_max)

    Lambda = args.Lambda
    if Lambda is None:
        Lambda = lambda_kinematics(args.p, args.xi, k_nom, th_nom, args.omega)
    domega = resonance_check(args.p, args.xi, k_nom, th_nom, args.omega,
                             args.n, args.B0)

    D_est = estimate_dxixi(args.W, args.D_diel, args.psi2, Lambda,
                           args.vg, args.vpar_costheta, args.dk)

    print("=" * 70)
    print("D_xixi ESTIMATE (analytic, single-mode delta approximation)")
    print("=" * 70)
    print(f"  phase-space point: p = {args.p}, xi = {args.xi}, n = {args.n}")
    print(f"  k = {k_nom:.2f} m^-1, theta_k = {th_nom:.3f} rad")
    print(f"  resonance check: Delta_omega = {domega:.3e} rad/s"
          + ("   (on resonance)" if abs(domega) < 0.05 * args.omega else "   (WARNING: off resonance)"))
    print(f"  Lambda = -xi + k_par v/omega = {Lambda:.3f}   (Lambda^2 = {Lambda**2:.2f})")
    print(f"  W = {args.W:.3e} J/m^3,  D_diel = {args.D_diel},  |psi|^2 = {args.psi2}")
    print(f"  |vg - v_par cos th| = {abs(args.vg - args.vpar_costheta):.3e} m/s,  "
          f"dk = {args.dk:.2f} m^-1")
    K = wave_space_factor(args.k_min, args.k_max, args.ktheta_min, args.ktheta_max)
    print(f"  k-space factor 2*pi*(kmax^3-kmin^3)/3*(cos th_min - cos th_max) = {K:.4e} m^-3")
    print(f"  -> amplitude matching W: A = W/(m c^2 K) = "
          f"{args.W / (M_E_C2_J * K):.4e}")
    print(f"\n  D_xixi_estimate = {D_est:.4e} s^-1")

    prog = load_program_dxixi(args.h5, args.p, args.xi, args.grid)
    if prog is None:
        print("\nNo program data to compare.")
        return

    print("\n" + "=" * 70)
    print("DREAM PROGRAM VALUE (net operator, QL + collisions)")
    print("=" * 70)
    print(f"  file: {args.h5}")
    print(f"  point used: p = {prog['p_actual']:.3f} (cell center), "
          f"xi = {prog['xi_actual']:.4f}")
    print(f"  Dxx at last saved time (t = {prog['t'][-1]:.4f} s): {prog['value']:.4e} s^-1")
    print(f"  max|Dxx| over the whole grid: {prog['max_all']:.4e} s^-1")

    print("\n" + "=" * 70)
    print("COMPARISON")
    print("=" * 70)
    if prog['value'] != 0:
        print(f"  D_est / D_prog = {D_est / prog['value']:.3e}")
    print(f"  D_est  = {D_est:.4e} s^-1")
    print(f"  D_prog = {prog['value']:.4e} s^-1")


if __name__ == '__main__':
    main()
