"""
可视化 DREAM 准线性扩散系数在相空间中的分布

使用方法:
    在 scripts 目录下运行：python plot_ql_diffusion_coefficients.py [选项]

命令行参数:
    --data_dir      数据文件目录 (默认: ../outputs)
    --plot_dir      图片保存目录 (默认: ../figures)
    --output_file   输出文件名前缀 (默认: ql_diffusion)

示例:
    python plot_ql_diffusion_coefficients.py
    python plot_ql_diffusion_coefficients.py --data_dir ../outputs --plot_dir ../figures_w
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import os
import sys
import argparse

# 添加 DREAM Python 库路径
sys.path.insert(0, '/data/zhzhou/DREAM/py')
from DREAM import DREAMOutput

def load_ql_diffusion_coefficients(data_dir, output_file='quasilinear_whistler_output.h5'):
    """
    从 DREAM 输出文件中加载准线性扩散系数
    
    Returns:
        D11, D22, D12: 扩散系数矩阵 (nr, np1_flux, np2) 或类似形状
        p_grid, xi_grid: 动量和投掷角网格
    """
    filepath = os.path.join(data_dir, output_file)
    
    if not os.path.exists(filepath):
        print(f"错误: 找不到文件 {filepath}")
        return None, None, None, None, None
    
    print(f"正在加载文件: {filepath}")
    
    try:
        do = DREAMOutput.DREAMOutput(filepath)
        
        # 获取最后一次时间步的数据
        last_timestep = len(do.time) - 1
        
        # 尝试从输出中读取扩散系数
        # 注意：DREAM 可能不会直接保存扩散系数，需要从方程项中提取
        # 这里我们尝试读取 f_re 的网格信息
        
        grid = do.getGrid('f_re')
        if grid is None:
            print("警告: 无法获取 f_re 网格信息")
            return None, None, None, None, None
        
        # 获取动量网格和投掷角网格
        p = grid.getP1()  # 动量网格中心
        xi = grid.getP2()  # 投掷角网格中心
        
        nr = grid.getNr()
        np1 = grid.getNp1()
        np2 = grid.getNp2()
        
        print(f"网格信息: nr={nr}, np1={np1}, np2={np2}")
        
        # 注意：DREAM 目前可能不直接保存扩散系数到输出文件
        # 我们需要通过其他方式获取，或者修改 C++ 代码来保存
        
        # 暂时返回网格信息，扩散系数需要特殊处理
        return None, None, None, p, xi
        
    except Exception as e:
        print(f"加载文件时出错: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None, None, None


def calculate_ql_coefficients_from_settings(settings_file, data_dir):
    """
    从设置文件重新计算准线性扩散系数（简化版本）
    
    这个方法通过读取波谱参数和等离子体参数，重新计算扩散系数
    """
    print("注意: 此功能需要访问 DREAM 内部的扩散系数计算")
    print("建议使用 C++ 代码中添加的诊断输出来获取准确的扩散系数")
    return None, None, None, None, None


def visualize_diffusion_coefficients(D11, D22, D12, p, xi, plot_dir, output_prefix='ql_diffusion'):
    """
    可视化扩散系数在相空间中的分布
    """
    if D11 is None or D22 is None or D12 is None:
        print("警告: 扩散系数数据为空，无法绘制")
        print("提示: 需要在 C++ 代码中添加扩散系数的诊断输出")
        return
    
    # 创建图形
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # 子图 1: D_pp (D11)
    ax1 = axes[0, 0]
    # 使用对数刻度，处理零值和负值
    D11_plot = np.abs(D11)
    D11_plot[D11_plot == 0] = np.min(D11_plot[D11_plot > 0]) if np.any(D11_plot > 0) else 1e-30
    
    im1 = ax1.pcolormesh(p, xi, D11_plot.T, cmap='viridis', norm=LogNorm(vmin=np.min(D11_plot), vmax=np.max(D11_plot)))
    ax1.set_xlabel('p (momentum)', fontsize=12)
    ax1.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
    ax1.set_title(r'$D_{pp}$ Distribution', fontsize=14, fontweight='bold')
    plt.colorbar(im1, ax=ax1, label=r'$D_{pp}$')
    ax1.grid(True, alpha=0.3)
    
    # 子图 2: D_xixi (D22)
    ax2 = axes[0, 1]
    D22_plot = np.abs(D22)
    D22_plot[D22_plot == 0] = np.min(D22_plot[D22_plot > 0]) if np.any(D22_plot > 0) else 1e-30
    
    im2 = ax2.pcolormesh(p, xi, D22_plot.T, cmap='plasma', norm=LogNorm(vmin=np.min(D22_plot), vmax=np.max(D22_plot)))
    ax2.set_xlabel('p (momentum)', fontsize=12)
    ax2.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
    ax2.set_title(r'$D_{\xi\xi}$ Distribution', fontsize=14, fontweight='bold')
    plt.colorbar(im2, ax=ax2, label=r'$D_{\xi\xi}$')
    ax2.grid(True, alpha=0.3)
    
    # 子图 3: D_pxi (D12)
    ax3 = axes[1, 0]
    D12_plot = np.abs(D12)
    D12_plot[D12_plot == 0] = np.min(D12_plot[D12_plot > 0]) if np.any(D12_plot > 0) else 1e-30
    
    im3 = ax3.pcolormesh(p, xi, D12_plot.T, cmap='inferno', norm=LogNorm(vmin=np.min(D12_plot), vmax=np.max(D12_plot)))
    ax3.set_xlabel('p (momentum)', fontsize=12)
    ax3.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
    ax3.set_title(r'$D_{p\xi}$ Distribution', fontsize=14, fontweight='bold')
    plt.colorbar(im3, ax=ax3, label=r'$D_{p\xi}$')
    ax3.grid(True, alpha=0.3)
    
    # 子图 4: 扩散系数比值 D_pxi / sqrt(D_pp * D_xixi)
    ax4 = axes[1, 1]
    # 这个比值应该在 [-1, 1] 范围内（柯西-施瓦茨不等式）
    denominator = np.sqrt(np.abs(D11 * D22))
    denominator[denominator == 0] = 1e-30
    ratio = D12 / denominator
    
    im4 = ax4.pcolormesh(p, xi, ratio.T, cmap='RdBu_r', vmin=-1, vmax=1)
    ax4.set_xlabel('p (momentum)', fontsize=12)
    ax4.set_ylabel(r'$\xi$ (pitch angle cosine)', fontsize=12)
    ax4.set_title(r'$D_{p\xi} / \sqrt{D_{pp} D_{\xi\xi}}$', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im4, ax=ax4)
    cbar.set_label('Correlation Coefficient')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 保存图片
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    output_path = os.path.join(plot_dir, f'{output_prefix}_coefficients.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"图片已保存: {output_path}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='可视化 DREAM 准线性扩散系数')
    parser.add_argument('--data_dir', type=str, default='../outputs', help='数据文件目录')
    parser.add_argument('--plot_dir', type=str, default='../figures', help='图片保存目录')
    parser.add_argument('--output_file', type=str, default='quasilinear_whistler_output.h5', help='输出文件名')
    parser.add_argument('--output_prefix', type=str, default='ql_diffusion', help='输出文件前缀')
    
    args = parser.parse_args()
    
    print("="*70)
    print("DREAM 准线性扩散系数可视化")
    print("="*70)
    
    # 尝试加载扩散系数
    D11, D22, D12, p, xi = load_ql_diffusion_coefficients(args.data_dir, args.output_file)
    
    if D11 is None:
        print("\n" + "="*70)
        print("重要提示:")
        print("="*70)
        print("当前 DREAM 输出文件不包含扩散系数的直接存储。")
        print("要获取扩散系数数据，需要：")
        print("1. 在 C++ 代码中添加诊断输出（已在 QuasilinearDiffusionTerm::Rebuild 中添加）")
        print("2. 或者使用 saveDiagnosticData() 方法保存扩散系数")
        print("3. 或者从运行时输出中提取验证信息")
        print("="*70)
        
        # 显示一个示例图（使用模拟数据）
        print("\n生成示例图（使用模拟数据）...")
        
        # 创建示例网格
        p_example = np.linspace(0, 50, 100)
        xi_example = np.linspace(-0.99, 0.99, 80)
        P, XI = np.meshgrid(p_example, xi_example)
        
        # 生成模拟的扩散系数（仅用于演示）
        D11_example = np.exp(-(P-20)**2/50) * np.exp(-(XI-0.3)**2/0.1) * 1e-2
        D22_example = np.exp(-(P-20)**2/50) * np.exp(-(XI-0.3)**2/0.1) * 1e-1
        D12_example = np.exp(-(P-20)**2/50) * np.exp(-(XI-0.3)**2/0.1) * 1e-3
        
        visualize_diffusion_coefficients(D11_example, D22_example, D12_example, 
                                        p_example, xi_example, args.plot_dir, args.output_prefix)
    else:
        # 使用真实数据绘制
        visualize_diffusion_coefficients(D11, D22, D12, p, xi, args.plot_dir, args.output_prefix)
    
    print("\n完成！")


if __name__ == "__main__":
    main()
