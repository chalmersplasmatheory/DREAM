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

# 获取分布函数数据 — 优先使用 f_re (动量范围更大，可达 p≫1)
if hasattr(do.eqsys, 'f_re') and do.eqsys.f_re is not None:
    f = do.eqsys.f_re[:, 0, :, :]   # [nt, nxi, np]
    p = do.grid.runaway.p[:]
    xi = do.grid.runaway.xi[:]
    print("Using f_re (runaway grid, extended momentum range)")
else:
    f = do.eqsys.f_hot[:, 0, :, :]  # [nt, nxi, np]
    p = do.grid.hottail.p[:]
    xi = do.grid.hottail.xi[:]
    print("Using f_hot (hottail grid)")

def _auto_ylim(values, pad_bottom=3, pad_top=3):
    """从数据中自动计算合适的 ylim，忽略 -inf。"""
    valid = np.isfinite(values)
    if not np.any(valid):
        return -20, 5
    vmin = np.min(values[valid])
    vmax = np.max(values[valid])
    return vmin - pad_bottom, vmax + pad_top


# ========= 图1: 所有时间步的 xi-积分分布函数叠加 =========
plt.figure(figsize=(6, 6))

# 选取指定时间点对应的索引
target_times = np.array([0.5, 0.9, 1.1, 1.2, 3.0])
selected_times = np.unique([np.argmin(np.abs(do.grid.t[:] - tt)) for tt in target_times])

colors = plt.cm.viridis(np.linspace(0, 1, len(selected_times)))
all_log = []

for i, t_idx in enumerate(selected_times):
    f_at_time = f[t_idx, :, :]
    f_vs_p = np.trapz(np.abs(f_at_time), xi, axis=0)
    f_nonzero = np.where(f_vs_p > 0, f_vs_p, 1e-30)
    log_f = np.log10(f_nonzero)
    all_log.append(log_f)
    
    plt.plot(p, log_f, color=colors[i], linewidth=2, 
             label=f't={do.grid.t[t_idx]:.3f} s')

plt.xlabel('p ($m_e c$)')
plt.ylabel(r'$\log_{10}( \int |f| \, d\xi )$')
plt.xlim(0, 25)
yl1, yh1 = _auto_ylim(np.concatenate(all_log))
plt.ylim(0, yh1)
plt.title(r'Distribution function evolution ($\xi$-integrated)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

file_path = os.path.join(args.plot_dir, 'f_time_evolution.png')
plt.savefig(file_path, dpi=300, bbox_inches='tight')
plt.show()

# ========= 图2: 4个时间点的详细分布 =========
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

num_times = f.shape[0]
detail_times = np.linspace(1, num_times-1, 4, dtype=int)  # t=0 时刻 f=0，跳过

for i, t_idx in enumerate(detail_times):
    ax = axes[i]
    f_at_time = f[t_idx, :, :]
    f_vs_p = np.trapz(np.abs(f_at_time), xi, axis=0)
    f_nonzero = np.where(f_vs_p > 0, f_vs_p, 1e-30)
    log_f = np.log10(f_nonzero)
    
    ax.plot(p, log_f, 'b-', linewidth=2)
    ax.set_xlabel('p ($m_e c$)')
    ax.set_ylabel(r'$\log_{10}(\int |f| \, d\xi)$')
    ax.set_xlim(0, 50)
    yl, yh = _auto_ylim(log_f)
    ax.set_ylim(yl, yh)
    ax.set_title(f't = {do.grid.t[t_idx]*1e6:.1f} μs')
    ax.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
file_path_detail = os.path.join(args.plot_dir, 'f_time_detail.png')
plt.savefig(file_path_detail, dpi=300, bbox_inches='tight')
plt.show()

# ========= 图3: 不同 pitch 角对比 =========
plt.figure(figsize=(6, 6))
nxi = xi.size
xi_indices = np.linspace(0, nxi-1, 5, dtype=int)
colors = plt.cm.plasma(np.linspace(0, 1, len(xi_indices)))
time_index = -1

all_log = []
for i, idx in enumerate(xi_indices):
    f_at_xi = f[time_index, idx, :]
    f_abs = np.abs(f_at_xi)
    f_nonzero = np.where(f_abs > 0, f_abs, 1e-30)
    log_f = np.log10(f_nonzero)
    all_log.append(log_f)
    
    plt.plot(p, log_f, color=colors[i], linewidth=2,
             label=f'ξ={xi[idx]:.2f}')

plt.xlabel('p ($m_e c$)')
plt.ylabel(r'$\log_{10}(|f|)$')
plt.xlim(0, 50)
yl, yh = _auto_ylim(np.concatenate(all_log))
plt.ylim(yl, yh)
plt.title(f'Distribution function at different pitch angles (t={do.grid.t[time_index]:.3f} s)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

file_path_xi = os.path.join(args.plot_dir, 'f_pitch_comparison.png')
plt.savefig(file_path_xi, dpi=300, bbox_inches='tight')
plt.show()

# ========= 图4: 捕获电子数量随能量的分布 =========
if hasattr(do.grid, 'xi0TrappedBoundary') and do.grid.xi0TrappedBoundary is not None:
    xi0 = do.grid.xi0TrappedBoundary[0]
    mask_trapped = np.abs(xi) < xi0

    # 能量 E_kin [MeV] = (sqrt(1+p²) - 1) * 0.511
    energy = (np.sqrt(1 + p**2) - 1) * 0.511

    plt.figure(figsize=(6, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_times)))
    all_log = []

    for i, t_idx in enumerate(selected_times):
        f_at_time = f[t_idx, :, :]          # [nxi, np]
        # 只对捕获区域积分
        f_trapped = np.trapz(np.abs(f_at_time[mask_trapped, :]), xi[mask_trapped], axis=0)
        f_nonzero = np.where(f_trapped > 0, f_trapped, 1e-30)
        log_f = np.log10(f_nonzero)
        all_log.append(log_f)

        plt.plot(energy, log_f, color=colors[i], linewidth=2,
                 label=f't={do.grid.t[t_idx]:.3f} s')

    plt.xlabel(r'Kinetic energy $E_k$ (MeV)')
    plt.ylabel(r'$\log_{10}( \int_{\rm trapped} |f| \, d\xi )$')
    plt.xlim(0, None)
    yl, yh = _auto_ylim(np.concatenate(all_log))
    plt.ylim(yl, yh)
    plt.title(r'Trapped electron distribution ($\xi$-integrated, $|\xi| < \xi_0$)')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)

    file_path_trapped = os.path.join(args.plot_dir, 'f_trapped_energy.png')
    plt.savefig(file_path_trapped, dpi=300, bbox_inches='tight')
    plt.show()

print(f"Plots saved to {args.plot_dir} directory")

# 关闭文件
do.close()