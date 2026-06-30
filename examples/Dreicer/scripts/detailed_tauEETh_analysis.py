#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')
from DREAM import DREAMOutput
from DREAM.Formulas.PlasmaParameters import getTauEETh

def main():
    # 加载DREAM输出
    do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')
    
    # 获取等离子体参数
    T = do.eqsys.T_cold[0, 0]  # 电子温度 (eV)
    n = do.eqsys.n_cold[0, 0]  # 电子密度 (m^-3)
    
    print(f"Plasma parameters from simulation:")
    print(f"  Temperature: {T:.2f} eV")
    print(f"  Density:     {n:.2e} m^-3")
    
    # 提取时间和tauEETh
    t = do.grid.t[1:]  # 时间网格（跳过第一个点）
    tauEETh_sim = do.other.fluid.tauEETh[:, 0]  # 热电子碰撞时间（模拟值）
    
    # 计算理论值
    tauEETh_theory = getTauEETh(T, n)
    
    print(f"Theoretical tauEETh: {tauEETh_theory:.3e} s")
    print(f"Average simulated tauEETh: {np.mean(tauEETh_sim):.3e} s")
    print(f"Difference: {abs(np.mean(tauEETh_sim) - tauEETh_theory)/tauEETh_theory*100:.2f}%")
    
    # 绘制tauEETh随时间的变化
    plt.figure(figsize=(10, 6))
    plt.plot(t * 1e6, tauEETh_sim, 'b-', linewidth=2, label='Simulated')
    plt.axhline(y=tauEETh_theory, color='r', linestyle='--', linewidth=2, 
                label=f'Theoretical = {tauEETh_theory:.3e} s')
    plt.axhline(y=np.mean(tauEETh_sim), color='g', linestyle=':', linewidth=2,
                label=f'Simulated average = {np.mean(tauEETh_sim):.3e} s')
    
    plt.xlabel('Time (μs)')
    plt.ylabel(r'Thermal electron collision time $\tau_{EE}^{Th}$ (s)')
    plt.title(r'Thermal electron collision time $\tau_{EE}^{Th}$ vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.yscale('log')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('../figures/detailed_tauEETh_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算并打印一些统计信息
    print(f"\nDetailed statistical summary for tauEETh:")
    print(f"  Mean:     {np.mean(tauEETh_sim):.3e} s")
    print(f"  Std dev:  {np.std(tauEETh_sim):.3e} s")
    print(f"  Min:      {np.min(tauEETh_sim):.3e} s")
    print(f"  Max:      {np.max(tauEETh_sim):.3e} s")
    print(f"  Relative variation: {np.std(tauEETh_sim)/np.mean(tauEETh_sim)*100:.2f}%")

if __name__ == "__main__":
    main()