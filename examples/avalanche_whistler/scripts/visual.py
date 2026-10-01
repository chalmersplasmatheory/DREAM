#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

# 添加DREAM路径
sys.path.append('/data/zhzhou/DREAM/py')

from DREAM import DREAMOutput

import argparse
import os

def load_simulation_results(data_dir=None):
    """
    加载所有模拟结果
    """
    if data_dir:
        # If a specific file is provided, compute E/Ec and use as key
        results = {}
        try:
            do = DREAMOutput(data_dir)
            # 从数据中提取 E/Ec
            E_field = do.eqsys.E_field.get()[0, 0]
            Ec_tot  = do.other.fluid.Ectot[:][0, 0]
            E_over_Ec = float(E_field / Ec_tot)
            results[E_over_Ec] = do
            print(f"Loaded data from {data_dir}, E/Ec = {E_over_Ec:.4f}")
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
    从模拟数据中计算增长率（物理量纲 + 导师归一化）
    
    导师归一化: gamma_norm = (nr/nrlast - 1)/dt * tauEERel * 2 * lnLambdaC
    
    Returns
    -------
    t_avg      : 1D array — 时间中点
    gamma      : 1D array — 物理增长率 (s^-1)
    gamma_norm : 1D array — 归一化增长率（无量纲，可对比 Rosenbluth）
    """
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]
    
    # 读取导师归一化所需的量
    tau = do.other.fluid.tauEERel[:][:, 0]
    lnL = do.other.fluid.lnLambdaC[:][:, 0]
    
    dt = np.diff(t)
    dn_re = np.diff(n_re)
    n_re_avg = 0.5 * (n_re[:-1] + n_re[1:])
    # tau/lnL 形状 (nt_other,) = (3000,)，已经对应时间中点
    # gamma 形状也是 (nt-1,) = (3000,)，直接对齐
    
    mask = n_re_avg > 1e5
    
    if np.sum(mask) > 0:
        gamma = np.where(mask, dn_re/dt/n_re_avg, 0)
        gamma_norm = np.where(mask, gamma * tau * 2.0 * lnL, 0)
        t_avg = 0.5 * (t[:-1] + t[1:])
        return t_avg[mask], gamma[mask], gamma_norm[mask]
    else:
        return np.array([]), np.array([]), np.array([])

def plot_growth_rates_vs_time(results, plot_dir):
    """
    绘制增长率随时间的变化
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for E_fac, do in results.items():
        t, gamma, _ = calculate_growth_rates(do)
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

def plot_normalized_growth_rate_vs_time(results, plot_dir):
    """
    绘制归一化增长率随时间的变化 + Rosenbluth 理论线
    归一化: gamma_norm = gamma_phys * tau * 2 * lnLambda
    Rosenbluth: gamma_norm_RP = 2 * (E/Ec - 1)
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    for ef, do in results.items():
        t, _, gamma_norm = calculate_growth_rates(do)
        if len(t) > 0:
            # 只绘制 t >= 0.5 的数据，避免早期异常值压扁纵轴
            mask_t = t >= 0.5
            ax.plot(t[mask_t], gamma_norm[mask_t], 'b-', linewidth=2, label=f'DREAM (E/Ec={ef:.2f})')
            
            # Rosenbluth 理论值: gamma_norm = 2*(E/Ec - 1)
            gamma_RP = 1.0 * (ef - 1.0)
            ax.axhline(y=gamma_RP, color='r', linestyle='--', linewidth=2,
                       label=f'Rosenbluth: $E/E_c-1$ = {gamma_RP:.2f}')
            print(f'  Rosenbluth prediction: 2*(E/Ec-1) = 2*({ef:.4f}-1) = {gamma_RP:.4f}')
    
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r"Normalized Growth Rate $\gamma \tau \cdot 2\ln\Lambda$")
    ax.set_title('Normalized Growth Rate vs Time')
    ax.legend(fontsize=12)
    ax.grid(True)
    ax.set_xlim(left=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_normalized_growth_rate_vs_time.png'), dpi=300)
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
    计算每个E/Ec值的平均增长率（物理 + 归一化）
    """
    E_factors = []
    avg_growth_rates = []
    avg_norm_rates = []
    
    for E_fac, do in results.items():
        t, gamma, gamma_norm = calculate_growth_rates(do)
        if len(gamma) > 10:
            avg_gamma = np.mean(gamma[-20:])
            avg_norm = np.mean(gamma_norm[-20:])
            E_factors.append(E_fac)
            avg_growth_rates.append(avg_gamma)
            avg_norm_rates.append(avg_norm)
            print(f"E/Ec = {E_fac}: gamma = {avg_gamma:.4e} s^-1  |  gamma_norm = {avg_norm:.4f}")
    
    return np.array(E_factors), np.array(avg_growth_rates), np.array(avg_norm_rates)
def plot_growth_rate_comparison(E_factors, avg_growth_rates, avg_norm_rates, plot_dir, do=None):
    """
    绘制平均增长率与E/Ec的关系（双图：物理量纲 + 归一化对比Rosenbluth）
    """
    # ---- 图1: 物理量纲 ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(E_factors, avg_growth_rates, 'bo-', label='DREAM Simulation')
    ax.set_xlabel('E/Ec')
    ax.set_ylabel('Average Growth Rate (s^-1)')
    ax.set_title('Average Dreicer Growth Rate vs Electric Field')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_avg_growth_rates.png'), dpi=300)
    plt.show()
    
    # ---- 图2: 归一化增长率 vs Rosenbluth ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    ax.plot(E_factors, avg_norm_rates, 'bo-', markersize=8,
            label='DREAM (normalized)')
    
    # Rosenbluth: gamma_norm = 2 * (E/Ec - 1)
    E_range = np.linspace(min(E_factors), max(E_factors), 100)
    rosenbluth = 2.0 * (E_range - 1.0)
    ax.plot(E_range, rosenbluth, 'r--', linewidth=2,
            label=r'Rosenbluth: $E/E_c-1$')
    
    ax.set_xlabel('E/Ec')
    ax.set_ylabel(r"Normalized Growth Rate $\gamma \tau \cdot 2\ln\Lambda$")
    ax.set_title('Normalized Growth Rate vs Rosenbluth Theory')
    ax.legend(fontsize=12)
    ax.grid(True)
    
    for i, (ef, nr) in enumerate(zip(E_factors, avg_norm_rates)):
        ax.annotate(f'{nr:.2f}', (ef, nr), textcoords='offset points',
                    xytext=(0, 10), ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'dreicer_normalized_vs_rosenbluth.png'), dpi=300)
    print(f"Normalized: last gamma_norm = {avg_norm_rates[-1]:.4f}")
    print(f"Rosenbluth: E/Ec-1     = {2*(E_factors[-1]-1):.4f}")
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
    
    # 绘制归一化增长率 vs 时间（含 Rosenbluth 理论线）
    print("Plotting normalized growth rate vs time...")
    plot_normalized_growth_rate_vs_time(results, args.plot_dir)
    
    # 计算并绘制平均增长率
    print("Calculating and plotting average growth rates...")
    E_factors, avg_growth_rates, avg_norm_rates = calculate_average_growth_rates(results)
    first_do = list(results.values())[0] if results else None
    plot_growth_rate_comparison(E_factors, avg_growth_rates, avg_norm_rates, args.plot_dir, do=first_do)
    
    print("Visualization completed. Plots saved as PNG files.")

if __name__ == "__main__":
    main()