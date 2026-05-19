#!/usr/bin/env python3
#
# Plot distribution function evolution at different times
#

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

sys.path.append('../../py')

from DREAM import *
from DREAM.DREAMOutput import DREAMOutput

import argparse

parser = argparse.ArgumentParser(description='Plot distribution function evolution')
parser.add_argument('--data_dir', type=str, default='../outputs/dreicer_with_fre_output.h5',
                    help='Path to the DREAM output HDF5 file')
parser.add_argument('--plot_dir', type=str, default='../figures',
                    help='Directory to save plots')
args = parser.parse_args()

# 创建保存图像的文件夹
os.makedirs(args.plot_dir, exist_ok=True)

# 加载数据
do = DREAMOutput(args.data_dir)

# 获取分布函数数据
f_hot = do.eqsys.f_hot[:, 0, :, :]  # 所有时间步，第一个径向点，所有xi和p
p = do.grid.hottail.p[:]
xi = do.grid.hottail.xi[:]

# 选择特定的pitch角索引（接近平行方向）
xi_index = -1  # 最后一个xi点，接近xi=1（平行方向）

# 创建一个图形，将所有时间步的线条画在同一张图中
plt.figure(figsize=(10, 6))

# 选择一些有代表性的时刻来绘制，避免线条过多看不清
# 这里我们每隔一定间隔选取一个时间步
num_times = f_hot.shape[0]
step_interval = max(1, num_times // 8)  # 最多绘制8条线
selected_times = range(0, num_times, step_interval)

# 为每个选定的时间步绘制分布函数
colors = plt.cm.viridis(np.linspace(0, 1, len(selected_times)))

for i, t_idx in enumerate(selected_times):
    f_at_time = f_hot[t_idx, :, :]  # 选定时间步，所有xi和p
    f_selected_xi = f_at_time[xi_index, :]  # 选定pitch角的分布函数
    # 避免取log(0)
    f_abs = np.abs(f_selected_xi)
    f_nonzero = np.where(f_abs > 0, f_abs, 1e-30)
    
    plt.plot(p, np.log10(f_nonzero), color=colors[i], linewidth=2, 
             label=f't={do.grid.t[t_idx]*1e6:.1f} μs')

plt.xlabel('p ($m_e c$)')
plt.ylabel(r'$\log_{10}(|f|)$')
plt.ylim(-20, 5)
plt.title('Distribution function time evolution at ξ≈1')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

# 保存图像到文件夹
file_path = os.path.join(args.plot_dir, 'f_time_evolution.png')
plt.savefig(file_path, dpi=300, bbox_inches='tight')
plt.show()

# 同时也创建一个单独时间步的高分辨率图
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

# 选择4个时间点进行详细展示
detail_times = np.linspace(0, num_times-1, 4, dtype=int)

for i, t_idx in enumerate(detail_times):
    ax = axes[i]
    f_at_time = f_hot[t_idx, :, :]
    f_selected_xi = f_at_time[xi_index, :]
    f_abs = np.abs(f_selected_xi)
    f_nonzero = np.where(f_abs > 0, f_abs, 1e-30)
    
    ax.plot(p, np.log10(f_nonzero), 'b-', linewidth=2)
    ax.set_xlabel('p ($m_e c$)')
    ax.set_ylabel(r'$\log_{10}(|f|)$')
    ax.set_ylim(-20, 5)
    ax.set_title(f't = {do.grid.t[t_idx]*1e6:.1f} μs')
    ax.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
file_path_detail = os.path.join(args.plot_dir, 'f_time_detail.png')
plt.savefig(file_path_detail, dpi=300, bbox_inches='tight')
plt.show()

# 创建一个展示不同pitch角度的对比图
plt.figure(figsize=(10, 6))
# 选择几个不同的xi值
xi_indices = [0, 5, 10, 15, 21]  # 从负到正的不同pitch角
colors = plt.cm.plasma(np.linspace(0, 1, len(xi_indices)))
time_index = -1  # 最后一个时间点

for i, idx in enumerate(xi_indices):
    f_at_xi = f_hot[time_index, idx, :]
    f_abs = np.abs(f_at_xi)
    f_nonzero = np.where(f_abs > 0, f_abs, 1e-30)
    
    plt.plot(p, np.log10(f_nonzero), color=colors[i], linewidth=2,
             label=f'ξ={xi[idx]:.2f}')

plt.xlabel('p ($m_e c$)')
plt.ylabel(r'$\log_{10}(|f|)$')
plt.ylim(-20, 5)
plt.title(f'Distribution function at different pitch angles (t={do.grid.t[time_index]*1e6:.1f} μs)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

file_path_xi = os.path.join(args.plot_dir, 'f_pitch_comparison.png')
plt.savefig(file_path_xi, dpi=300, bbox_inches='tight')
plt.show()

print(f"Plots saved to {args.plot_dir} directory")

# 关闭文件
do.close()