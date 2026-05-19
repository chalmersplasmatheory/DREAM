# DREAM Kinetic Equation Implementation

## 完整的动理学方程

DREAM求解的bounce-averaged动理学方程为：

$$\frac{\partial f}{\partial t} + \nabla \cdot (\mathbf{A} f) - \nabla \cdot (\mathsf{D} \nabla f) = S$$

其中 $\mathbf{A}$ 是advection系数，$\mathsf{D}$ 是diffusion系数，$S$ 是源项。

---

## 1. 电场加速项 (Electric Field Term)

### **源码位置**
- 实现：`/data/zhzhou/DREAM/src/Equations/Kinetic/ElectricFieldTerm.cpp`
- 头文件：`/data/zhzhou/DREAM/include/DREAM/Equations/Kinetic/ElectricFieldTerm.hpp`

### **物理方程**

原始Vlasov方程中的电场项：
$$eE_\parallel \xi \frac{\partial f}{\partial p} + \frac{eE_\parallel}{p}(1-\xi^2)\frac{\partial f}{\partial \xi}$$

Bounce average后（P-XI网格）：

**动量方向 advection** ($F_1$):
```cpp
// ElectricFieldTerm.cpp line 59
F1(ir, i, j) += xi0 * E_xi_bounceAvg_f1;
```

对应方程：
$$A^p = \xi_0 \cdot \frac{e}{m_e c} \frac{\{E_\parallel \xi\}}{\xi_0} \cdot \frac{\sqrt{\langle B^2 \rangle}}{\langle B \rangle}$$

**Pitch角方向 advection** ($F_2$):
```cpp
// ElectricFieldTerm.cpp line 72
F2(ir, i, j) += E_xi_bounceAvg_f2 * (1-xi0*xi0)/p;
```

对应方程：
$$A^{\xi_0} = \frac{e}{m_e c} \frac{\{E_\parallel \xi\}}{\xi_0} \cdot \frac{\sqrt{\langle B^2 \rangle}}{\langle B \rangle} \cdot \frac{1-\xi_0^2}{p}$$

### **关键计算**

```cpp
// Line 49-55: 计算归一化因子
sqrtB2OverB = sqrt(grid->GetRadialGrid()->GetFSA_B2(ir)) / grid->GetRadialGrid()->GetFSA_B(ir);
E = Constants::ec * E_term[ir] /(Constants::me * Constants::c);

// Line 55: 计算 {Eξ}/ξ₀ × √⟨B²⟩/⟨B⟩
E_xi_bounceAvg_f1 = E * xiAvgTerm_f1[j*(np1+1)+i] * sqrtB2OverB;
```

其中：
- `xiAvgTerm_f1` = $\{E_\parallel \xi\}/\xi_0$ （通过bounce average预计算）
- `sqrtB2OverB` = $\sqrt{\langle B^2 \rangle}/\langle B \rangle$ （几何因子）

### **Bounce Average函数定义**

```cpp
// RadialGrid.hpp line 58-59
static real_t BA_FUNC_XI(real_t xiOverXi0, real_t, real_t, real_t, void*)
    {return xiOverXi0;}  // X(θ) = ξ/ξ₀

static constexpr int_t BA_PARAM_XI[5] = {1,0,0,0,1};
// [POWER_XI, POWER_B, POWER_R, POWER_NABLA_R, normalization]
```

---

## 2. 同步辐射阻尼 (Synchrotron Radiation)

### **源码位置**
- 实现：`/data/zhzhou/DREAM/src/Equations/Kinetic/SynchrotronTerm.cpp`
- 理论推导：`/data/zhzhou/DREAM/doc/notes/theory.tex` (line 1891-1923)

### **物理方程**

同步辐射引起的动量空间advection：

**原始形式**（theory.tex line 1893-1896）：
$$\mathbf{A}_s = -\frac{1}{\gamma \tau_s}\left(p_\perp\frac{\partial \mathbf{p}}{\partial p_\perp} + p\frac{p_\perp^2}{m_e^2c^2}\frac{\partial \mathbf{p}}{\partial p}\right)$$

其中特征时间：
$$\frac{1}{\tau_s} = \frac{e^4 B^2}{6\pi\varepsilon_0 m_e^3c^3}$$

**在 (p, ξ₀) 坐标下**（theory.tex line 1914-1919）：

动量方向：
$$A_s^p = -\frac{p\gamma}{\tau_{s,\min}} (1-\xi_0^2)\frac{B^3}{B_{\min}^3}$$

Pitch角方向：
$$A_s^{\xi_0} = \frac{1-\xi_0^2}{\gamma\tau_{s,\min}} \frac{B^2}{B_{\min}^2}\frac{\xi^2}{\xi_0}$$

展开为两项：
$$A_s^{\xi_0} = \frac{\xi_0(1-\xi_0^2)}{\gamma \tau_{s,\min}} \frac{B^3}{B_{\min}^3} - \frac{1}{\gamma \tau_{s,\min}}\frac{1-\xi_0^2}{\xi_0}\frac{B^2}{B_{\min}^2}\left(\frac{B}{B_{\min}}-1\right)$$

### **代码实现**

```cpp
// SynchrotronTerm.cpp line 103: P方向 advection
return -preFactor * p*sqrt(1+p*p)*(1-xi0*xi0) * BA1_f1[j*(mg->GetNp1()+1)+i];
```

对应：
$$F_1 = -\underbrace{B_{\min}^2 \cdot \text{constPrefactor}}_{\text{preFactor}} \cdot p\gamma(1-\xi_0^2) \cdot \underbrace{\left\{\frac{B^3}{B_{\min}^3}\right\}}_{\text{BA1\_f1}}$$

```cpp
// SynchrotronTerm.cpp line 113: ξ₀方向 advection
return +preFactor * (1-xi0*xi0)*xi0/gamma * BA2_f2[j*mg->GetNp1()+i];
```

对应：
$$F_2 = +B_{\min}^2 \cdot \text{constPrefactor} \cdot \frac{(1-\xi_0^2)\xi_0}{\gamma} \cdot \underbrace{\left\{\frac{B^2}{B_{\min}^2}\frac{\xi^2}{\xi_0^2}\right\}}_{\text{BA2\_f2}}$$

### **Bounce Average函数定义**

```cpp
// RadialGrid.hpp line 62-65
static real_t BA_FUNC_B_CUBED(real_t, real_t BOverBmin, real_t, real_t, void*)
    {return BOverBmin*BOverBmin*BOverBmin;}  // X = (B/B_min)³

static real_t BA_FUNC_XI_SQUARED_B_SQUARED(real_t xiOverXi0, real_t BOverBmin, real_t, real_t, void*)
    {return BOverBmin*BOverBmin*xiOverXi0*xiOverXi0;}  // X = (B/B_min)²·(ξ/ξ₀)²

static constexpr int_t 
    BA_PARAM_B_CUBED[5] = {0,3,0,0,1},           // [0, 3, 0, 0, 1]
    BA_PARAM_XI_SQUARED_B_SQUARED[5] = {2,2,0,0,1};  // [2, 2, 0, 0, 1]
```

---

## 3. Fokker-Planck 碰撞算符

FP碰撞算符包含三部分：**慢化 (Slowing Down)**、**能量扩散 (Energy Diffusion)**、**Pitch-angle散射 (Pitch Scattering)**

### **3.1 Pitch-angle散射项**

#### **源码位置**
- 实现：`/data/zhzhou/DREAM/src/Equations/Kinetic/PitchScatterTerm.cpp`

#### **物理方程**

原始FP碰撞算符的pitch散射部分：
$$C_{\text{pitch}} = \frac{\nu_D}{2} \frac{\partial}{\partial \xi}\left[(1-\xi^2)\frac{\partial f}{\partial \xi}\right]$$

其中 $\nu_D$ 是deflection频率。

**Bounce average后**（P-XI网格）：

Diffusion系数 $D_{22}$：
```cpp
// PitchScatterTerm.cpp line 61
D22(ir,i,j) += commonFactor_f2 * (1 - xi0*xi0);
```

对应方程：
$$D^{\xi_0\xi_0} = \frac{1}{2} \nu_D \cdot \underbrace{\left\{\frac{B_{\min}}{B}\frac{\xi^2}{\xi_0^2}\right\}}_{\text{xiBAvg\_f2}} \cdot (1-\xi_0^2)$$

#### **Bounce Average函数**

```cpp
// RadialGrid.hpp line 60-61
static real_t BA_FUNC_XI_SQUARED_OVER_B(real_t xiOverXi0, real_t BOverBmin, real_t, real_t, void*)
    {return xiOverXi0*xiOverXi0/BOverBmin;}  // X = (ξ/ξ₀)²/(B/B_min)

static constexpr int_t BA_PARAM_XI_SQUARED_OVER_B[5] = {2,-1,0,0,1};
// [POWER_XI=2, POWER_B=-1, POWER_R=0, POWER_NABLA_R=0, norm=1]
```

---

### **3.2 慢化项 (Slowing Down / Friction)**

#### **源码位置**
- 实现：`/data/zhzhou/DREAM/src/Equations/Kinetic/SlowingDownTerm.cpp`

#### **物理方程**

原始FP碰撞算符的摩擦部分：
$$C_{\text{friction}} = -\frac{1}{p^2}\frac{\partial}{\partial p}\left(p^2 \nu_s p f\right)$$

其中 $\nu_s$ 是slowing-down频率。

**Bounce average后**（与角度无关，直接advection）：

动量方向：
```cpp
// SlowingDownTerm.cpp line 46
F1(ir, i, j) -= mg->GetP1_f(i) * nu_s_f1[ir][j*(n1[ir]+1)+i];
```

对应方程：
$$A^p = -p \cdot \nu_s$$

**注意**：慢化项不依赖角度，因此不需要bounce average！

---

### **3.3 能量扩散项 (Energy Diffusion)**

#### **源码位置**
- 实现：`/data/zhzhou/DREAM/src/Equations/Kinetic/EnergyDiffusionTerm.cpp`

#### **物理方程**

原始FP碰撞算符的能量扩散部分：
$$C_{\text{energy}} = \frac{1}{p^2}\frac{\partial}{\partial p}\left(p^2 D^{pp} \frac{\partial f}{\partial p}\right)$$

其中 $D^{pp} = m_e T_{\text{cold}} \gamma \nu_s$（theory.tex line 1940-1942）。

**Bounce average后**：

Diffusion系数：
$$D^{pp} = (m_e c)^2 \nu_\parallel = (m_e c)^2 \frac{T_{\text{cold}}}{m_e c^2} \gamma \nu_s$$

同样**不依赖角度**，无需bounce average。

---

## 4. Bounce Average参数化系统

### **通用函数签名**

```cpp
real_t F(real_t xiOverXi0, real_t BOverBmin, real_t ROverR0, real_t NablaR2, void* par)
```

### **幂次参数表示**

DREAM使用整数数组指定幂次组合：
```cpp
int_t Flist[5] = {POWER_XI, POWER_B, POWER_R, POWER_NABLA_R, normalization_flag}
```

### **预定义的Bounce Average函数**

| 函数名 | 表达式 | 幂次参数 | 用途 |
|--------|--------|---------|------|
| `BA_FUNC_UNITY` | 1 | `{0,0,0,0,1}` | 相空间体积 Vp |
| `BA_FUNC_XI` | ξ/ξ₀ | `{1,0,0,0,1}` | 电场加速 {Eξ} |
| `BA_FUNC_XI_SQUARED_OVER_B` | (ξ/ξ₀)²/(B/B_min) | `{2,-1,0,0,1}` | Pitch散射 {ν_D ξ²/B} |
| `BA_FUNC_B_CUBED` | (B/B_min)³ | `{0,3,0,0,1}` | 同步辐射 {B³} |
| `BA_FUNC_XI_SQUARED_B_SQUARED` | (B/B_min)²·(ξ/ξ₀)² | `{2,2,0,0,1}` | 同步辐射 {B²ξ²/ξ₀²} |

### **如何添加新的物理项**

1. **定义Bounce Average函数**：
```cpp
static real_t BA_FUNC_MY_PHYSICS(real_t xiOvXi0, real_t BOvBmin, real_t ROvR0, real_t NablaR2, void*)
    {return pow(xiOvXi0, a) * pow(BOvBmin, b) * pow(ROvR0, c) * pow(NablaR2, d);}
```

2. **定义幂次参数**：
```cpp
static constexpr int_t BA_PARAM_MY_PHYSICS[5] = {a, b, c, d, 1};
```

3. **在Grid.cpp中预计算**：
```cpp
SetBounceAveragePXi(BA_my_physics_f1, FLUXGRIDTYPE_P1, 
                    RadialGrid::BA_FUNC_MY_PHYSICS, nullptr, 
                    RadialGrid::BA_PARAM_MY_PHYSICS);
```

4. **在物理项中使用**：
```cpp
const real_t *BA_my = grid->GetBA_my_physics_f1(ir);
F1(ir, i, j) += coefficient * BA_my[index];
```

---

## 5. 完整方程总结

对于P-XI网格，DREAM求解的bounce-averaged动理学方程为：

$$\frac{\partial f}{\partial t} + \frac{1}{p^2}\frac{\partial}{\partial p}\left[p^2 (A^p f - D^{pp}\frac{\partial f}{\partial p})\right] + \frac{1}{p}\frac{\partial}{\partial \xi_0}\left[(1-\xi_0^2)(A^{\xi_0} f - D^{\xi_0\xi_0}\frac{\partial f}{\partial \xi_0})\right] = S$$

其中各项系数：

| 系数 | 表达式 | 是否bounce averaged |
|------|--------|-------------------|
| $A^p_{E}$ | $\xi_0 \frac{e}{m_ec} \frac{\{E\xi\}}{\xi_0} \frac{\sqrt{\langle B^2\rangle}}{\langle B\rangle}$ | ✅ |
| $A^{\xi_0}_{E}$ | $\frac{e}{m_ec} \frac{\{E\xi\}}{\xi_0} \frac{\sqrt{\langle B^2\rangle}}{\langle B\rangle} \frac{1-\xi_0^2}{p}$ | ✅ |
| $A^p_{s}$ | $-\frac{p\gamma}{\tau_{s,\min}}(1-\xi_0^2)\{\frac{B^3}{B_{\min}^3}\}$ | ✅ |
| $A^{\xi_0}_{s}$ | $\frac{(1-\xi_0^2)\xi_0}{\gamma\tau_{s,\min}}\{\frac{B^2}{B_{\min}^2}\frac{\xi^2}{\xi_0^2}\}$ | ✅ |
| $A^p_{\text{fric}}$ | $-p\nu_s$ | ❌ (无角度依赖) |
| $D^{pp}$ | $(m_ec)^2 \frac{T}{m_ec^2}\gamma\nu_s$ | ❌ (无角度依赖) |
| $D^{\xi_0\xi_0}$ | $\frac{1}{2}\nu_D \{\frac{B_{\min}}{B}\frac{\xi^2}{\xi_0^2}\}(1-\xi_0^2)$ | ✅ |

---

## 6. 关键洞察

### **为什么有些项需要bounce average，有些不需要？**

- **需要bounce average**：系数依赖于极向角θ（通过B(θ), ξ(θ), R(θ)等）
  - 电场加速：E_∥可能随θ变化，ξ(θ)显式依赖θ
  - 同步辐射：强烈依赖B(θ)
  - Pitch散射：ν_D和ξ都依赖θ

- **不需要bounce average**：系数与θ无关或已解析积分
  - 慢化项：ν_s只依赖p和等离子体参数，与θ无关
  - 能量扩散：D^{pp}只依赖p和T_cold，与θ无关

### **Bounce Average的核心价值**

通过将快变的极向运动平均掉，DREAM可以：
1. **大幅减少计算量**：无需分辨极向时间尺度
2. **保持物理精度**：精确捕捉捕获/通行粒子的差异
3. **灵活扩展**：通过函数指针机制轻松添加新物理

这就是DREAM方程构造的完整图景！
