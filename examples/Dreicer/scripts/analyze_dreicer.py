#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

# 添加DREAM路径
sys.path.append('/data/zhzhou/DREAM/py/')

from DREAM import DREAMOutput

def analyze_dreicer_generation():
    """
    分析Dreicer生成率
    """
    # 加载模拟结果
    do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')
    
    # 提取时间网格和物理量
    t = do.grid.t[:]
    n_re = do.eqsys.n_re.get()[:, 0]
    
    # 计算Dreicer生成率 (dn_re/dt)
    dt = np.diff(t)
    dn_re = np.diff(n_re)
    
    # 避免除零，只在n_re足够大时计算增长率
    n_re_avg = 0.5 * (n_re[1:] + n_re[:-1])
    mask = n_re_avg > 1e5  # 只考虑n_re > 1e5的情况
    
    gamma_dreicer = np.zeros_like(dn_re)
    gamma_dreicer[mask] = dn_re[mask]/dt[mask]/n_re_avg[mask]
    
    # 时间平均值
    t_avg = 0.5 * (t[1:] + t[:-1])
    
    # 只绘制0.2微秒之后的数据
    time_mask = t_avg >= 0.2e-6  # 0.2微秒 = 0.2e-6秒
    combined_mask = mask & time_mask
    
    # 绘制逃逸电子密度随时间的变化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1.plot(t, n_re, 'b-', linewidth=2)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Runaway electron density (m$^{-3}$)')
    ax1.set_title('Runaway Electron Density vs Time')
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.set_yscale('log')
    
    # 绘制Dreicer增长率（仅0.2微秒之后的数据）
    ax2.plot(t_avg[combined_mask], gamma_dreicer[combined_mask], 'r-', linewidth=2)
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Dreicer Growth Rate (s$^{-1}$)')
    ax2.set_title('Dreicer Growth Rate vs Time (t >= 0.2 μs)')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('../figures/dreicer_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算平均增长率（仅使用有效的数据点）
    valid_gammas = gamma_dreicer[combined_mask]
    if len(valid_gammas) > 10:  # 确保有足够的数据点
        # 去除前10%和后10%以避免瞬态效应
        start_idx = len(valid_gammas) // 10
        end_idx = 9 * len(valid_gammas) // 10
        avg_growth_rate = np.mean(valid_gammas[start_idx:end_idx])
        std_growth_rate = np.std(valid_gammas[start_idx:end_idx])
    else:
        avg_growth_rate = np.mean(valid_gammas)
        std_growth_rate = np.std(valid_gammas)
    
    print(f"Average Dreicer growth rate: {avg_growth_rate:.3e} s^-1")
    print(f"Standard deviation: {std_growth_rate:.3e} s^-1")
    
    # 与理论值比较（如果需要的话）
    E = 6  # V/m
    Ec = 30  # V/m (大约值，对于n=5e19 m^-3, T=100 eV)
    E_ratio = E / Ec
    
    print(f"Electric field ratio E/Ec: {E_ratio:.2f}")
    
    return avg_growth_rate, E_ratio

def main():
    """
    主函数
    """
    print("Analyzing Dreicer generation results...")
    
    try:
        avg_rate, E_ratio = analyze_dreicer_generation()
        print("Analysis completed. Results saved to '../figures/dreicer_analysis.png'")
    except Exception as e:
        print(f"Error during analysis: {e}")
        raise

if __name__ == "__main__":
    main()