#!/usr/bin/env python3
#
# Batch processor for DREAM-SOFT Green function multiplication
# Reads f_re from DREAM HDF5 output, interpolates to SOFT grid,
# multiplies with Green function, and saves results per timestep.
#
'''
python dream_green_processor.py <dream_output.h5> <green.mat> <output_dir> [--radial-index 0] [--max-steps N]

Example:
cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts && python dream_green_processor.py \
    ../outputs/quasilinear_whistler_output.h5 \
    /path/to/green.mat \
    ../outputs/green_results \
    --radial-index 0
'''

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import glob
import sys
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import argparse

sys.path.append('/data/zhzhou/DREAM/py/')


def tos(v):
    """Convert HDF5 string to Python string"""
    return "".join(map(chr, v[:,:][:,0].tolist()))


def _normalize_thetap(thetap, green_func, axis):
    """Normalize thetap axis to xi=cos(thetap) ascending order.

    Returns (XI, green_func) where XI is strictly ascending in [-1, 1].
    Handles: degree/radian detection, range validation, arbitrary input order.
    """
    thetap = np.array(thetap, dtype=float)
    gf = np.array(green_func, dtype=float)

    if thetap.max() > 2 * np.pi + 1e-6:
        print(f"  Warning: thetap max={thetap.max():.4f} > 2pi, assuming DEGREES")
        thetap = np.deg2rad(thetap)

    if thetap.min() < -1e-9 or thetap.max() > np.pi + 1e-9:
        raise ValueError(
            f"thetap range [{thetap.min():.4f}, {thetap.max():.4f}] outside [0, pi]; "
            "xi=cos(thetap) is not injective, cannot convert safely.")

    order = np.argsort(thetap)
    thetap = thetap[order]
    gf = np.take(gf, order, axis=axis)

    XI = np.cos(thetap)
    reverse_idx = np.arange(len(XI))[::-1]
    XI = XI[::-1]
    gf = np.take(gf, reverse_idx, axis=axis)

    assert np.all(np.diff(XI) > 0), "XI not ascending after normalization"
    return XI, gf


def load_green_function(filename):
    """Load Green's function from SOFT output"""
    print(f"Loading Green's function from: {filename}")
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Green's function file not found: {filename}")

    with h5py.File(filename, 'r') as f:
        green_func = f['func'][:]
        par1 = f['param1'][:]
        par2 = f['param2'][:]

        frmt = tos(f['type'])
        par1n = tos(f['param1name'])
        par2n = tos(f['param2name'])

        if frmt != '12':
            raise Exception("Invalid format of Green's function: '{0}'.".format(frmt))

        P, XI = None, None

        if (par1n == 'p' and par2n == 'xi'):
            P, XI = par1, par2
        elif (par1n == 'xi' and par2n == 'p'):
            XI, P = par1, par2
            green_func = green_func.T
        elif (par1n == 'p' and par2n == 'thetap'):
            P = par1
            XI, green_func = _normalize_thetap(par2, green_func, axis=1)
            print(f"  Converted thetap -> xi: range [{par2.min():.4f}, {par2.max():.4f}] -> [{XI.min():.4f}, {XI.max():.4f}]")
        elif (par1n == 'thetap' and par2n == 'p'):
            P = par2
            XI, green_func = _normalize_thetap(par1, green_func, axis=0)
            green_func = green_func.T
            print(f"  Converted thetap -> xi: range [{par1.min():.4f}, {par1.max():.4f}] -> [{XI.min():.4f}, {XI.max():.4f}]")
        else:
            raise Exception("Expected p-xi or p-thetap parameters, found: '{0}' & '{1}'.".format(par1n, par2n))

        if green_func.shape != (P.size, XI.size):
            if green_func.shape == (XI.size, P.size):
                print(f"  Transposing green_func {green_func.shape} -> {(P.size, XI.size)}")
                green_func = green_func.T
            else:
                raise ValueError(
                    f"green_func shape {green_func.shape} inconsistent with "
                    f"len(P)={P.size}, len(XI)={XI.size}")

    print(f"Green function shape: {green_func.shape}")
    return P, XI, green_func


def load_dream_data(dream_h5_path, radial_index=0):
    """Load f_re distribution from DREAM HDF5 output

    Parameters:
    dream_h5_path: Path to DREAM output HDF5 file
    radial_index: Which radial grid point to use (default 0)

    Returns:
    time_array: numpy 1D array of time values (Nt,)
    dream_p: numpy 1D array of momentum grid (Np,)
    dream_xi: numpy 1D array of pitch angle grid (Nxi,)
    f_all: numpy 3D array of distribution (Nt, Nxi, Np)
    """
    from DREAM.DREAMOutput import DREAMOutput

    print(f"Loading DREAM data from: {dream_h5_path}")
    print(f"  Using radial index: {radial_index}")

    do = DREAMOutput(dream_h5_path, lazy=True)

    time_array = do.grid.t[:]
    dream_p = do.grid.runaway.p[:]
    dream_xi = do.grid.runaway.xi[:]
    f_all = do.eqsys.f_re[:, radial_index, :, :]  # (Nt, Nxi, Np)

    do.close()

    print(f"  Loaded {len(time_array)} time steps")
    print(f"  p grid: {len(dream_p)} points, range [{dream_p.min():.4f}, {dream_p.max():.4f}]")
    print(f"  xi grid: {len(dream_xi)} points, range [{dream_xi.min():.4f}, {dream_xi.max():.4f}]")
    print(f"  f_re shape per step: ({len(dream_xi)}, {len(dream_p)})")

    return time_array, dream_p, dream_xi, f_all


def interpolate_dream_to_soft(dream_p, dream_xi, dream_f_2d, soft_p_grid, soft_xi_grid):
    """Interpolate DREAM distribution to SOFT grid using structured interpolation

    Parameters:
    dream_p: 1D array of DREAM momentum grid (Np_dream,)
    dream_xi: 1D array of DREAM xi grid (Nxi_dream,)
    dream_f_2d: 2D array of distribution (Nxi_dream, Np_dream) — DREAM native layout
    soft_p_grid: 1D array of SOFT momentum grid (Np_soft,)
    soft_xi_grid: 1D array of SOFT xi grid (Nxi_soft,)

    Returns:
    2D array (Np_soft, Nxi_soft) matching Green function convention
    """
    # Build interpolator on DREAM's (xi, p) grid
    interpolator = RegularGridInterpolator(
        (dream_xi, dream_p), dream_f_2d,
        method='linear', bounds_error=False, fill_value=0.0
    )

    # Create SOFT evaluation grid — meshgrid with indexing='ij' gives (Np, Nxi)
    soft_p_mesh, soft_xi_mesh = np.meshgrid(soft_p_grid, soft_xi_grid, indexing='ij')

    # Evaluate: need (Nxi, Np) points for the interpolator
    # soft_xi_mesh has shape (Np, Nxi), transpose to get (Nxi, Np) pairs
    eval_points = np.column_stack([soft_xi_mesh.ravel(), soft_p_mesh.ravel()])
    result = interpolator(eval_points).reshape(soft_p_mesh.shape)

    return result  # (Np_soft, Nxi_soft)


def multiply_and_save(green_func, distribution, green_p, green_xi, time_val, output_dir, timestep_idx, no_plot=False):
    """Multiply Green function with distribution and save result with plot"""

    if green_func.shape != distribution.shape:
        raise ValueError(f"Shape mismatch: green {green_func.shape} vs f {distribution.shape}. "
                         "Check interpolation output and Green function loading.")

    result = green_func * distribution

    data_dir = os.path.join(output_dir, "data")
    plots_dir = os.path.join(output_dir, "plots")
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    Path(plots_dir).mkdir(parents=True, exist_ok=True)

    output_file = os.path.join(data_dir, f"green_times_f_timestep_{timestep_idx:04d}.h5")
    with h5py.File(output_file, 'w') as f:
        f['p'] = green_p
        f['xi'] = green_xi
        f['result'] = result
        f['time'] = time_val
        f.attrs['description'] = 'Green function multiplied with distribution function'
        f.attrs['format_version'] = '1.0'
        f.attrs['green_function_convention'] = 'power_per_dp_dxi'

    plot_file = None
    if not no_plot:
        plot_file = _generate_plot(green_func, distribution, green_p, green_xi, result, time_val, plots_dir, timestep_idx)

    return output_file, plot_file


def _generate_plot(green_func, distribution, green_p, green_xi, result, time_val, plots_dir, timestep_idx):
    """Generate 3-panel plot for a single timestep."""
    plot_file = os.path.join(plots_dir, f"green_times_f_t_{time_val:.4f}.png")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Panel 1: Heatmap of Green*f in (p, xi)
    im1 = axes[0].pcolormesh(green_p, green_xi, result.T, shading='auto', cmap='viridis')
    axes[0].set_xlabel('Momentum p')
    axes[0].set_ylabel(r'Pitch angle $\xi = \cos\theta$')
    axes[0].set_title(f'Green x f at t={time_val:.4f}')
    plt.colorbar(im1, ax=axes[0], label='Value')

    # Panel 2: Integrated over pitch angle
    integrand = result * green_p[:, np.newaxis]**2
    integrated_over_xi = np.trapz(integrand, x=green_xi, axis=1)
    axes[1].semilogy(green_p, integrated_over_xi, 'b-', linewidth=2)
    axes[1].set_xlabel('Momentum p')
    axes[1].set_ylabel(r'$\int G \cdot f \cdot p^2 \, d\xi$')
    axes[1].set_title(f'Integrated at t={time_val:.4f}')
    axes[1].grid(True, alpha=0.3)

    # Panel 3: p_par-p_perp weight map
    P, XI = np.meshgrid(green_p, green_xi, indexing='ij')
    PPAR = P * XI
    PPERP = P * np.sqrt(np.clip(1 - XI**2, 0, None))
    F_weighted = result * PPERP
    fmax = np.amax(np.abs(F_weighted))
    if fmax > 0:
        F_weighted = F_weighted / fmax
    axes[2].contourf(PPAR, PPERP, F_weighted, levels=50, cmap='hot')
    axes[2].set_xlabel(r'$p_\parallel / mc$')
    axes[2].set_ylabel(r'$p_\perp / mc$')
    axes[2].set_title(f'Green x f weighted (p_par, p_perp) at t={time_val:.4f}')
    axes[2].set_aspect('equal')

    plt.tight_layout()
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return plot_file


def clean_previous_results(output_dir):
    """Remove per-timestep files left by a previous run in the same directory.

    The output file index (0000, 0001, ...) restarts on every run, so leftover
    files from a run with a different time selection are not overwritten and
    would later be mixed into the power calculation.
    """
    removed = 0
    patterns = [
        (os.path.join(output_dir, "data"), "green_times_f_timestep_*.h5"),
        (os.path.join(output_dir, "plots"), "green_times_f_t_*.png"),
    ]
    for directory, pattern in patterns:
        for filename in glob.glob(os.path.join(directory, pattern)):
            os.remove(filename)
            removed += 1

    if removed:
        print(f"Removed {removed} file(s) from a previous run in: {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="DREAM-SOFT Green function batch processor")
    parser.add_argument("dream_h5_path", help="Path to DREAM output HDF5 file")
    parser.add_argument("green_mat_path", help="Path to SOFT green.mat file")
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument("--radial-index", type=int, default=0, help="Radial grid index (default: 0)")
    parser.add_argument("--dt", type=float, default=None, help="Time interval for subsampling in seconds (default: use all steps)")
    parser.add_argument("--max-steps", type=int, default=None, help="Maximum number of time steps to process")
    parser.add_argument("--t-min", type=float, default=None, help="Minimum time to process")
    parser.add_argument("--t-max", type=float, default=None, help="Maximum time to process")
    parser.add_argument("--no-plot", action="store_true", help="Skip plot generation (save HDF5 only)")
    parser.add_argument("--no-clean", action="store_true",
                        help="Keep per-timestep files from a previous run in output_dir "
                             "(default: remove them)")

    args = parser.parse_args()

    print("=" * 60)
    print("DREAM-SOFT GREEN FUNCTION BATCH PROCESSOR")
    print("=" * 60)
    print(f"DREAM output: {args.dream_h5_path}")
    print(f"Green function: {args.green_mat_path}")
    print(f"Output directory: {args.output_dir}")
    print(f"Radial index: {args.radial_index}")
    print(f"Max steps: {args.max_steps if args.max_steps else 'All'}")
    print(f"Time range: [{args.t_min if args.t_min else '-inf'}, {args.t_max if args.t_max else 'inf'}]")
    print(f"dt: {args.dt if args.dt else 'all steps'}")
    print(f"Plot: {'disabled' if args.no_plot else 'enabled'}")
    print("=" * 60)

    # Load Green's function
    green_p, green_xi, green_func = load_green_function(args.green_mat_path)

    # Load DREAM data
    time_array, dream_p, dream_xi, f_all = load_dream_data(args.dream_h5_path, args.radial_index)

    # Filter time steps by range
    if args.t_min is not None or args.t_max is not None:
        mask = np.ones(len(time_array), dtype=bool)
        if args.t_min is not None:
            mask &= time_array >= args.t_min
        if args.t_max is not None:
            mask &= time_array <= args.t_max
        indices = np.where(mask)[0]
        time_array = time_array[indices]
        f_all = f_all[indices]
        print(f"Filtered to {len(time_array)} time steps within time range")

    # Subsample by dt
    if args.dt is not None:
        t_min_val = time_array[0]
        t_max_val = time_array[-1]
        target_times = np.arange(t_min_val, t_max_val + args.dt * 0.5, args.dt)
        indices = [np.argmin(np.abs(time_array - tt)) for tt in target_times]
        indices = np.unique(indices)
        time_array = time_array[indices]
        f_all = f_all[indices]
        print(f"Subsampled to {len(time_array)} time steps (dt={args.dt})")

    # Apply max-steps limit
    if args.max_steps is not None and len(time_array) > args.max_steps:
        time_array = time_array[:args.max_steps]
        f_all = f_all[:args.max_steps]
        print(f"Limited to first {args.max_steps} time steps")

    # Remove files from a previous run only after the inputs loaded, so a bad
    # path does not wipe earlier results.
    if not args.no_clean:
        clean_previous_results(args.output_dir)

    print(f"\nStarting batch processing of {len(time_array)} time steps...")
    successful_steps = 0

    for i, time_val in enumerate(time_array):
        try:
            print(f"Step {i+1}/{len(time_array)}: Processing time={time_val:.4f}")

            # Extract distribution for this timestep: (Nxi, Np) → transpose → (Np, Nxi)
            dream_f_2d = f_all[i, :, :]  # (Nxi, Np)

            # Interpolate from DREAM grid to SOFT grid
            soft_dist = interpolate_dream_to_soft(dream_p, dream_xi, dream_f_2d, green_p, green_xi)

            # Multiply and save
            data_file, plot_file = multiply_and_save(
                green_func, soft_dist, green_p, green_xi, time_val,
                args.output_dir, i, no_plot=args.no_plot
            )

            print(f"  Saved: {os.path.basename(data_file)}")
            if plot_file:
                print(f"  Plot:  {os.path.basename(plot_file)}")
            successful_steps += 1

        except Exception as e:
            print(f"  Error processing step {i}: {e}")
            continue

    print("\n" + "=" * 60)
    print(f"BATCH PROCESSING COMPLETED")
    print(f"Successfully processed: {successful_steps}/{len(time_array)} time steps")
    print(f"Results saved to: {args.output_dir}")
    print(f"Data files: {args.output_dir}/data/")
    print(f"Plot files: {args.output_dir}/plots/")
    print("=" * 60)


if __name__ == "__main__":
    main()
