#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append('/data/zhzhou/DREAM/py/')
from DREAM import DREAMOutput

def main():
    # 加载DREAM输出
    do = DREAMOutput('../outputs/dreicer_with_fre_output.h5')
    
    # 提取时间和tauEETh
    t = do.grid.t[1:]  # 时间网格（跳过第一个点，因为other.quantity不包含t=0的数据）
    tauEETh = do.other.fluid.tauEETh[:, 0]  # 热电子碰撞时间（去掉额外的维度）
    
    print(f"tauEETh shape: {tauEETh.shape}")
    print(f"Time grid shape: {t.shape}")
    print(f"First few tauEETh values:\n{tauEETh[:10]}")
    print(f"First few time values:\n{t[:10]}")
    
    # 绘制tauEETh随时间的变化
    plt.figure(figsize=(10, 6))
    plt.plot(t, tauEETh, 'b-', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel(r'Thermal electron collision time $\tau_{EE}^{Th}$ (s)')
    plt.title(r'Thermal electron collision time $\tau_{EE}^{Th}$ vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.yscale('log')
    
    # 添加文本注释显示一些统计信息
    avg_tauEETh = np.mean(tauEETh)
    plt.axhline(y=avg_tauEETh, color='r', linestyle='--', 
                label=f'Average = {avg_tauEETh:.2e} s')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('../figures/tauEETh_vs_time.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 计算并打印一些统计信息
    print(f"\nStatistical summary for tauEETh:")
    print(f"  Mean:     {np.mean(tauEETh):.3e} s")
    print(f"  Std dev:  {np.std(tauEETh):.3e} s")
    print(f"  Min:      {np.min(tauEETh):.3e} s")
    print(f"  Max:      {np.max(tauEETh):.3e} s")
    print(f"  Relative variation: {np.std(tauEETh)/np.mean(tauEETh)*100:.2f}%")

if __name__ == "__main__":
    main()