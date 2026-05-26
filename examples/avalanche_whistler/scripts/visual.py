#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

# 添加DREAM路径
sys.path.append('../../py')

from DREAM import DREAMOutput

import argparse
import os

def load_simulation_results(data_dir=None):
    """
    加载所有模拟结果
    """
    if data_dir:
        # If a specific file is provided, just load that one
        results = {}
        try:
            do = DREAMOutput(data_dir)
            results[data_dir] = do
            print(f"Loaded data from {data_dir}")
        except Exception as e:
            print(f"Failed to load data from {data_dir}: {e}")
        return results

    E_factors = [1.0, 1.5, 2.0, 3.0]  # 跳过E/Ec = 0.5的数据
    results = {}
    
    for E_fac in E_factors:
        filename = f'dreicer_validation_Efac{E_fac}.h5'
        if Path(filename).exists():
            try:
                do = DREAMOutput(filename)
                results[E_fac] = do
                print(f"Loaded data for E/Ec = {E_fac}")
            except Exception as e:
                print(f"Failed to load data for E/Ec = {E_fac}: {e}")
        else:
            print(f"File not found: {filename}")
            
    return results

def calculate_growth_rates(do):
    """
    从模拟数据中计算增长率
    """
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]  # 取第一个径向点的数据
    
    # 计算增长率: gamma = (1/n_re) * dn_re/dt
    # 使用对数导数: gamma = d(ln n_re)/dt
    dt = np.diff(t)
    dn_re = np.diff(n_re)
    
    # 避免除以零或非常小的数
    n_re_avg = 0.5 * (n_re[:-1] + n_re[1:])
    
    # 只有在n_re足够大时才计算增长率
    mask = n_re_avg > 1e5
    
    if np.sum(mask) > 0:
        gamma = np.where(mask, dn_re/dt/n_re_avg, 0)
        t_avg = 0.5 * (t[:-1] + t[1:])
        return t_avg[mask], gamma[mask]
    else:
        return np.array([]), np.array([])

def plot_growth_rates_vs_time(results, plot_dir):
    """
    绘制增长率随时间的变化
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for E_fac, do in results.items():
        t, gamma = calculate_growth_rates(do)
        if len(t) > 0:
            ax.plot(t, gamma, label=f'E/Ec = {E_fac}')
    
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Dreicer Growth Rate (s⁻¹)')
    ax.set_title('Dreicer Growth Rate vs Time')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_growth_rate_vs_time.png'), dpi=300)
    plt.show()

def plot_n_re_vs_time(results, plot_dir):
    """
    绘制逃逸电子密度随时间的变化
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for E_fac, do in results.items():
        t = do.grid.t[:]
        n_re = do.eqsys.n_re.get()[:, 0]  # 取第一个径向点的数据
        ax.plot(t, n_re, label=f'E/Ec = {E_fac}')
    
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Runaway Electron Density (m⁻³)')
    ax.set_title('Runaway Electron Density vs Time')
    ax.set_yscale('log')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_n_re_vs_time.png'), dpi=300)
    plt.show()

def calculate_average_growth_rates(results):
    """
    计算每个E/Ec值的平均增长率
    """
    E_factors = []
    avg_growth_rates = []
    
    for E_fac, do in results.items():
        t, gamma = calculate_growth_rates(do)
        if len(gamma) > 10:  # 确保有足够的数据点
            # 取最后20个时间步的平均值
            avg_gamma = np.mean(gamma[-20:])
            E_factors.append(E_fac)
            avg_growth_rates.append(avg_gamma)
            label_str = str(E_fac)
            print(f"E/Ec = {label_str}: Average growth rate = {avg_gamma:.4e} s⁻¹")
    
    return np.array(E_factors), np.array(avg_growth_rates)

def plot_growth_rate_comparison(E_factors, avg_growth_rates, plot_dir):
    """
    绘制平均增长率与E/Ec的关系
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    ax.plot(E_factors, avg_growth_rates, 'bo-', label='DREAM Simulation')
    ax.set_xlabel('E/Ec')
    ax.set_ylabel('Average Dreicer Growth Rate (s⁻¹)')
    ax.set_title('Average Dreicer Growth Rate vs Electric Field')
    ax.legend()
    ax.grid(True)
    
    # 添加第二个y轴显示归一化增长率
    ax2 = ax.twinx()
    if len(avg_growth_rates) > 0:
        normalized_rates = avg_growth_rates/np.max(np.abs(avg_growth_rates))
        ax2.plot(E_factors, normalized_rates, 'ro--', label='Normalized')
        ax2.set_ylabel('Normalized Growth Rate')
        ax2.legend(loc='lower right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_avg_growth_rates.png'), dpi=300)
    plt.show()

def main():
    """
    主函数
    """
    parser = argparse.ArgumentParser(description='Visualize DREAM simulation results')
    parser.add_argument('--data_dir', type=str, default=None,
                        help='Path to the DREAM output HDF5 file (optional)')
    parser.add_argument('--plot_dir', type=str, default='../figures',
                        help='Directory to save plots')
    args = parser.parse_args()

    os.makedirs(args.plot_dir, exist_ok=True)

    print("Loading DREAM simulation results...")
    results = load_simulation_results(args.data_dir)
    
    if not results:
        print("No simulation results found!")
        return
    
    print(f"Successfully loaded {len(results)} simulation results")
    
    # 绘制逃逸电子密度随时间的变化
    print("Plotting runaway electron density vs time...")
    plot_n_re_vs_time(results, args.plot_dir)
    
    # 绘制增长率随时间的变化
    print("Plotting growth rates vs time...")
    plot_growth_rates_vs_time(results, args.plot_dir)
    
    # 计算并绘制平均增长率
    print("Calculating and plotting average growth rates...")
    E_factors, avg_growth_rates = calculate_average_growth_rates(results)
    plot_growth_rate_comparison(E_factors, avg_growth_rates, args.plot_dir)
    
    print("Visualization completed. Plots saved as PNG files.")

if __name__ == "__main__":
    main()