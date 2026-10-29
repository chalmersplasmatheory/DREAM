#!/usr/bin/env python3
#
# Calculate synchrotron radiation total power from Green function multiplication results
# Compute total power by integrating Green×f over momentum space for each time step
#
'''
# Prerequisite: run dream_green_processor.py first to produce per-timestep radiation data

# Basic calculation
python dream_se_power.py /data/results/dir

# Set the plot time range
python dream_se_power.py /data/path --t-min 0.5 --t-max 2.0

# Mark wave injection intervals
python dream_se_power.py /data/path --wave-intervals "3.4,3.5;3.6,3.7;3.8,3.9"
'''

import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import glob
from pathlib import Path
import argparse


def calculate_total_power_from_hdf5(hdf5_file):
    """
    Calculate total synchrotron power from HDF5 file containing Green×f result

    Returns:
    tuple: (time_value, total_power)
    """
    with h5py.File(hdf5_file, 'r') as f:
        p_grid = f['p'][:]
        xi_grid = f['xi'][:]
        result = f['result'][:]
        time_val = f['time'][()]

        # Total power = 2π ∫∫ G·f·p² dp dξ
        integrand = result * p_grid[:, np.newaxis]**2
        integral_over_xi = np.trapz(integrand, x=xi_grid, axis=1)
        total_power = 2 * np.pi * np.trapz(integral_over_xi, x=p_grid)

        return time_val, total_power


def process_all_timesteps(results_dir):
    """
    Process all timestep files and calculate total power for each

    Returns:
    tuple: (times, powers) arrays
    """
    data_dir = os.path.join(results_dir, "data")
    hdf5_files = sorted(glob.glob(os.path.join(data_dir, "green_times_f_timestep_*.h5")))

    if not hdf5_files:
        raise FileNotFoundError(f"No HDF5 files found in {data_dir}")

    print(f"Found {len(hdf5_files)} timestep files")

    times = []
    powers = []

    for i, hdf5_file in enumerate(hdf5_files):
        try:
            time_val, total_power = calculate_total_power_from_hdf5(hdf5_file)
            times.append(time_val)
            powers.append(total_power)
            print(f"Step {i+1}/{len(hdf5_files)}: t={time_val:.4f}, Power={total_power:.3e}")
        except Exception as e:
            print(f"Error processing {hdf5_file}: {e}")
            continue

    return np.array(times), np.array(powers)


def parse_wave_intervals(intervals_str):
    """Parse wave interval string like '3.4,3.5;3.6,3.7' into list of (start, end) tuples"""
    if not intervals_str:
        return None
    intervals = []
    for pair in intervals_str.split(';'):
        pair = pair.strip()
        if not pair:
            continue
        parts = pair.split(',')
        if len(parts) != 2:
            raise ValueError(f"Invalid interval format: '{pair}'. Expected 'start,end'.")
        intervals.append((float(parts[0]), float(parts[1])))
    return intervals


def plot_power_vs_time(times, powers, output_dir, t_min=None, t_max=None, wave_intervals=None):
    """
    Plot total synchrotron power vs time and save the figure
    """
    mask = np.ones(len(times), dtype=bool)
    if t_min is not None:
        mask &= times >= t_min
    if t_max is not None:
        mask &= times <= t_max
    plot_times = times[mask]
    plot_powers = powers[mask]

    if len(plot_times) == 0:
        print("Warning: no data points in the specified time range!")
        return None

    plots_dir = os.path.join(output_dir, "plots")
    Path(plots_dir).mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(5, 5))

    plt.plot(plot_times, plot_powers, 'b-', linewidth=2, marker='o', markersize=4)

    if wave_intervals:
        actual_t_min = plot_times[0]
        actual_t_max = plot_times[-1]
        for i, (t_start, t_end) in enumerate(wave_intervals):
            if t_start <= actual_t_max and t_end >= actual_t_min:
                plt.axvspan(t_start, t_end, alpha=0.2, color='gray',
                          label='Wave Injection' if i == 0 else '')
        if wave_intervals:
            plt.legend(loc='best')

    plt.xlabel('Time t')
    plt.ylabel('Total Synchrotron Power')
    title = 'Synchrotron Radiation Total Power vs Time'
    if t_min is not None or t_max is not None:
        t_range = f"t = [{t_min if t_min is not None else ''}, {t_max if t_max is not None else ''}]"
        title += f'\n{t_range}'
    plt.title(title)
    plt.grid(True, alpha=0.3)

    plot_file = os.path.join(plots_dir, "synchrotron_total_power_vs_time.png")
    if wave_intervals:
        plot_file = plot_file.replace('.png', '_wave.png')
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Plot saved to: {plot_file}")
    return plot_file


def save_power_data(times, powers, output_dir):
    """Save power data to CSV file"""
    data_dir = os.path.join(output_dir, "data")
    Path(data_dir).mkdir(parents=True, exist_ok=True)

    csv_file = os.path.join(data_dir, "synchrotron_total_power.csv")
    with open(csv_file, 'w') as f:
        f.write("time,power\n")
        for t, p in zip(times, powers):
            f.write(f"{t:.6f},{p:.6e}\n")

    print(f"Data saved to: {csv_file}")
    return csv_file


def main():
    parser = argparse.ArgumentParser(description="Calculate synchrotron radiation total power from Green function multiplication results")
    parser.add_argument("results_dir", help="Directory containing HDF5 result files (with data/ subfolder)")
    parser.add_argument("--t-min", type=float, default=None, help="Minimum time for plot range")
    parser.add_argument("--t-max", type=float, default=None, help="Maximum time for plot range")
    parser.add_argument("--wave-intervals", type=str, default=None,
                        help='Wave injection intervals to shade, format: "3.4,3.5;3.6,3.7;3.8,3.9"')

    args = parser.parse_args()

    wave_intervals = parse_wave_intervals(args.wave_intervals)
    if wave_intervals:
        print(f"Wave injection intervals: {wave_intervals}")

    print("=" * 70)
    print("SYNCHROTRON RADIATION TOTAL POWER CALCULATION")
    print("=" * 70)
    print(f"Results directory: {args.results_dir}")
    print("=" * 70)

    try:
        times, powers = process_all_timesteps(args.results_dir)

        if len(times) == 0:
            print("No valid data found!")
            return

        print(f"\nPower Statistics:")
        print(f"  Time range: {times[0]:.4f} to {times[-1]:.4f}")
        print(f"  Power range: {np.min(powers):.3e} to {np.max(powers):.3e}")
        print(f"  Mean power: {np.mean(powers):.3e}")
        print(f"  Peak power at t={times[np.argmax(powers)]:.4f}")

        plot_file = plot_power_vs_time(times, powers, args.results_dir,
                                       t_min=args.t_min, t_max=args.t_max,
                                       wave_intervals=wave_intervals)

        csv_file = save_power_data(times, powers, args.results_dir)

        print("\n" + "=" * 70)
        print("CALCULATION COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print(f"Processed {len(times)} time steps")
        print(f"Plot saved: {plot_file}")
        print(f"Data saved: {csv_file}")
        print("=" * 70)

    except Exception as e:
        print(f"Error: {e}")
        return


if __name__ == "__main__":
    main()
