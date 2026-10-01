import numpy as np
import matplotlib.pyplot as plt

# ================== 物理参数 ==================
xi_prime = 0.999         # 初级电子俯仰角余弦（用于第一张图）
gamma_prime = 16.0       # 初级电子洛伦兹因子
pCutAvalanche = 2.0      # DREAM 中的动量截断
gamma_min = np.sqrt(1 + pCutAvalanche**2)  # γ_min = √(1+2²) = √5 ≈ 2.236

# 托卡马克环径比
epsilon = 0.24            # a/R (典型值)
# xi_trapped = sqrt(1 - Bmin/Bmax) = sqrt(2ε/(1+ε))
xi_trapped = np.sqrt(2 * epsilon / (1 + epsilon))

# 运动学允许的最大γ (能量较低电子定义为次级)
gamma_max = (gamma_prime + 1) / 2.0   # = 8.5
if gamma_max <= gamma_min:
    raise ValueError("No allowed secondary energy range")

# ================== 辅助函数 ==================
def p_from_gamma(gamma):
    """动量 p = √(γ² - 1)"""
    return np.sqrt(gamma**2 - 1)

def dsigma_dgamma(gamma_prime, gamma):
    """Møller微分截面 dσ/dγ  (相对值，忽略2πr₀²等常数)"""
    nu = (gamma - 1) / (gamma_prime - 1)
    if nu <= 0 or nu >= 1:
        return 0.0
    x = 1.0 / (nu * (1.0 - nu))
    prefactor = gamma_prime**2 / ((gamma_prime - 1)**3 * (gamma_prime + 1))
    bracket = x**2 - 3*x + ((gamma_prime - 1) / gamma_prime)**2 * (1 + x)
    return prefactor * bracket

def compute_Pi_params(gamma_prime, xi_prime, gamma):
    """计算ξ₁, ξ₂ (公式10b, 10c)"""
    xi1 = xi_prime * np.sqrt((gamma_prime + 1) * (gamma - 1) /
                              ((gamma_prime - 1) * (gamma + 1)))
    xi2 = np.sqrt(2 * (gamma_prime - gamma) /
                   ((gamma_prime - 1) * (gamma + 1)) * (1 - xi_prime**2))
    return xi1, xi2

def Pi_value(xi, xi1, xi2):
    """Π(ξ) 值：在区间内为1/(π√(ξ₂²-(ξ-ξ₁)²))，否则0"""
    if xi2 <= 0:
        return 1e6 if abs(xi - xi1) < 1e-8 else 0.0
    if abs(xi - xi1) <= xi2:
        return 1.0 / (np.pi * np.sqrt(xi2**2 - (xi - xi1)**2))
    else:
        return 0.0

def F_pass(gamma_prime, xi_prime, gamma, xi_trapped):
    """
    计算次级电子落在通行区 (|ξ| >= xi_trapped) 的概率分数。
    即对 Π(ξ) 在通行区内解析积分，归一化到总概率1。
    """
    xi1, xi2 = compute_Pi_params(gamma_prime, xi_prime, gamma)
    if xi2 <= 0:
        return 1.0 if abs(xi1) >= xi_trapped else 0.0
    low = max(xi1 - xi2, -1.0)
    high = min(xi1 + xi2, 1.0)

    def I_seg(a, b):
        """∫_a^b Π dξ"""
        if b <= a:
            return 0.0
        arg_a = max(-1.0, min(1.0, (a - xi1) / xi2))
        arg_b = max(-1.0, min(1.0, (b - xi1) / xi2))
        return (np.arcsin(arg_b) - np.arcsin(arg_a)) / np.pi

    total = 0.0
    # 右通行区: ξ >= xi_trapped
    if high >= xi_trapped:
        a = max(low, xi_trapped)
        b = high
        total += I_seg(a, b)
    # 左通行区: ξ <= -xi_trapped
    if low <= -xi_trapped:
        a = low
        b = min(high, -xi_trapped)
        total += I_seg(a, b)
    return max(min(total, 1.0), 0.0)

# ================== 构建网格 ==================
n_gamma = 200            # γ 方向网格数
n_xi = 300               # ξ 方向网格数

gamma_vals = np.linspace(gamma_min, gamma_max, n_gamma)
xi_vals = np.linspace(-1.0, 1.0, n_xi)

# 存储概率密度 (相对值)
density = np.zeros((n_gamma, n_xi))

for i, g in enumerate(gamma_vals):
    dsig = dsigma_dgamma(gamma_prime, g)
    if dsig <= 0:
        continue
    xi1, xi2 = compute_Pi_params(gamma_prime, xi_prime, g)
    for j, xi in enumerate(xi_vals):
        if xi2 == 0 and abs(xi - xi1) > 1e-8:
            continue
        prob = Pi_value(xi, xi1, xi2)
        weight = gamma_vals[i]**2 / (gamma_vals[i]**2 - 1) if g > 1 else 0
        density[i, j] = dsig * prob * weight

# ================== 统计 ==================
total_intensity = np.sum(density)
nonzero_pixels = np.sum(density > 0)
peak_val = np.max(density)
peak_idx = np.unravel_index(np.argmax(density), density.shape)
peak_gamma = gamma_vals[peak_idx[0]]
peak_xi = xi_vals[peak_idx[1]]

print("="*70)
print("Møller Scattering Analysis with DREAM pCutAvalanche")
print("="*70)
print(f"Primary electron:   γ' = {gamma_prime},  p' = {p_from_gamma(gamma_prime):.2f}")
print(f"pCutAvalanche:       p_cut = {pCutAvalanche}")
print(f"Minimum γ:           γ_min = √(1+p_cut²) = {gamma_min:.4f}")
print(f"Maximum γ:           γ_max = (γ'+1)/2 = {gamma_max:.3f}")
print(f"Integration range:   γ ∈ [{gamma_min:.3f}, {gamma_max:.3f}]")
print(f"                     p ∈ [{p_from_gamma(gamma_min):.3f}, {p_from_gamma(gamma_max):.3f}]")
print(f"Toroidal:            ε = a/R = {epsilon}")
print(f"Trapped boundary:    ξ_trapped = √(2ε/(1+ε)) = {xi_trapped:.4f}")
print("-"*70)
print(f"Grid: {n_gamma}×{n_xi} = {n_gamma*n_xi}")
print(f"Non-zero pixels: {nonzero_pixels} ({100*nonzero_pixels/(n_gamma*n_xi):.1f}%)")
print(f"Peak at: γ = {peak_gamma:.3f}, ξ = {peak_xi:.4f}")
print(f"Peak value: {peak_val:.4e} (arb.)")
print(f"ξ₁ at γ=γ_max: {compute_Pi_params(gamma_prime, xi_prime, gamma_max)[0]:.4f}")
print(f"ξ₂ at γ=γ_max: {compute_Pi_params(gamma_prime, xi_prime, gamma_max)[1]:.4f}")

# ================== 图1：次级电子二维分布 ==================
plt.figure(figsize=(7, 7))
vmax = np.percentile(density, 99.5)
vmin = 0

# 只取 ξ ≥ 0 的部分
n_xi_half = n_xi // 2
xi_half = xi_vals[n_xi_half:]
extent = [gamma_min, gamma_max, 0, 1]
plt.imshow(density[:, n_xi_half:].T, origin='lower', aspect='auto', extent=extent,
           cmap='inferno', vmin=vmin, vmax=vmax)
plt.colorbar(label='Relative production rate (arb. units)')
plt.xlabel(r'Secondary $\gamma$', fontsize=14)
plt.ylabel(r'Secondary $\xi = p_\parallel/p$', fontsize=14)
plt.title(rf'Secondary electron source from Møller scattering'
          '\n'
          rf'$\gamma\'={gamma_prime}$, $\xi\'={xi_prime}$, '
          rf'$\varepsilon={epsilon}$, $\gamma_{{\min}}={gamma_min:.2f}$', fontsize=14)

# ξ₁(γ) 曲线
xi1_curve = np.array([compute_Pi_params(gamma_prime, xi_prime, g)[0] for g in gamma_vals])
plt.plot(gamma_vals, xi1_curve, 'c--', linewidth=1.5, label=r'$\xi_1(\gamma)$')

# ±ξ₂ 包络线
xi2_curve = np.array([compute_Pi_params(gamma_prime, xi_prime, g)[1] for g in gamma_vals])
plt.plot(gamma_vals, xi1_curve + xi2_curve, 'w:', linewidth=1.0, alpha=0.7)
plt.plot(gamma_vals, xi1_curve - xi2_curve, 'w:', linewidth=1.0, alpha=0.7,
         label=r'$\xi_1 \pm \xi_2$')

# 捕获边界线
plt.axhline(xi_trapped, color='r', linestyle='-', linewidth=2, alpha=0.6)
plt.text((gamma_min + gamma_max) / 2, xi_trapped,
         rf'$\xi_{{\rm trapped}} = {xi_trapped:.3f}$',
         color='r', fontsize=13, fontweight='bold', ha='center', va='bottom',
         bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
# 捕获区标注
plt.gca().fill_between([gamma_min, gamma_max], 0, xi_trapped, 
                        color='red', alpha=0.06)
plt.text(gamma_min+0.05, xi_trapped/2, 'trapped\n(no avalanche)', 
         color='r', fontsize=11, fontweight='bold', alpha=0.6)
# 积分范围标注
plt.text(gamma_min+0.05, 0.92, 
         rf'Integration: $\gamma \in [{gamma_min:.3f}, {gamma_max:.3f}]$',
         fontsize=11, color='#374151',
         bbox=dict(facecolor='white', edgecolor='#D8E1EA', alpha=0.8))

plt.legend(fontsize=11, loc='lower left')
plt.tight_layout()
plt.savefig('MecDavidMoller_output.png', dpi=150)
print(f"\nSaved MecDavidMoller_output.png")

# ================== 图2：通行概率 vs 初级投掷角 ==================
print("\n========== 扫描初级电子投掷角 ξ' ==========")
print(f"Integration range: γ ∈ [{gamma_min:.3f}, {gamma_max:.3f}] (pCut={pCutAvalanche})")
print(f"xi_trapped = {xi_trapped:.4f} (ε = {epsilon})")

# 次级电子 γ 方向上的截面和权重
dsig_arr = np.array([dsigma_dgamma(gamma_prime, g) for g in gamma_vals])
weight_arr = gamma_vals**2 / (gamma_vals**2 - 1.0)

# 分母（通量归一化，与 ξ' 无关）
denom = np.trapz(dsig_arr * weight_arr, gamma_vals)

# 扫描 ξ'
xi_prime_scan = np.linspace(0.5, 0.9999, 300)
f_pass_scan = np.empty_like(xi_prime_scan)

for k, xip in enumerate(xi_prime_scan):
    fp_vals = np.array([F_pass(gamma_prime, xip, g, xi_trapped) for g in gamma_vals])
    numer = np.trapz(dsig_arr * weight_arr * fp_vals, gamma_vals)
    f_pass_scan[k] = numer / denom

print(f"ξ' scan range: [{xi_prime_scan[0]:.3f}, {xi_prime_scan[-1]:.5f}], {len(xi_prime_scan)} points")
print(f"F_pass range:  [{f_pass_scan.min()*100:.1f}%, {f_pass_scan.max()*100:.1f}%]")

print("\nKey data points:")
for target in [0.90, 0.95, 0.98, 0.99, 0.999, 0.9999]:
    idx = np.argmin(np.abs(xi_prime_scan - target))
    actual = xi_prime_scan[idx]
    print(f"  ξ'={target:.4f} (1-ξ'={1-actual:.2e}): F_pass = {f_pass_scan[idx]*100:.3f}%")

# 特殊询问的点：1-ξ'=0.05 (ξ'=0.95)
val_005 = np.interp(0.05, 1.0 - xi_prime_scan[::-1], f_pass_scan[::-1]) * 100
print(f"\n>> 1-ξ' = 0.05 (ξ'=0.95): F_pass = {val_005:.2f}% <<")

# ========== 图2：通行概率 vs 初级投掷角 ==========
plt.figure(figsize=(8, 7))

plt.plot(xi_prime_scan, f_pass_scan * 100, 'b-', lw=2.5, label=r'$F_{\rm pass}(\xi\')$')

# 标注关键 ξ' 值
for target in [0.6, 0.7, 0.8, 0.85, 0.9, 0.95]:
    idx = np.argmin(np.abs(xi_prime_scan - target))
    y = f_pass_scan[idx] * 100
    plt.axvline(target, color='gray', linestyle=':', alpha=0.3)
    plt.text(target+0.005, y+2, rf"$\xi'={target:.2g}$", fontsize=10, ha='left', fontweight='bold')

# RP 极限
plt.axhline(100, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
plt.text(0.995, 102, "RP limit ($\\xi'=1$): 100%", 
         color='gray', fontsize=11, fontweight='bold', ha='right')

# 标注积分范围
plt.text(0.52, 15, 
         rf'Integration: $\gamma \in [{gamma_min:.2f}, {gamma_max:.2f}]$',
         fontsize=11, color='#374151',
         bbox=dict(facecolor='#F8FAFC', edgecolor='#D8E1EA', alpha=0.8))

plt.xlabel(r"Primary $\xi' = \cos\theta'$", fontsize=14)
plt.ylabel(r"Passing fraction $F_{\rm pass}$ (%)", fontsize=14)
plt.title(rf"Avalanche efficiency vs seed pitch angle"
          rf"  ($\gamma\'={gamma_prime},\; \varepsilon={epsilon}$)", fontsize=14)
plt.ylim(0, 110)
plt.xlim(0.5, 1.0)
plt.grid(alpha=0.3)
plt.legend(fontsize=12, loc='lower left')
plt.tight_layout()
plt.savefig('MecDavidMoller_passing_fraction.png', dpi=150)
print(f"\nSaved MecDavidMoller_passing_fraction.png")
print("="*70)