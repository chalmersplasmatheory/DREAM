# DREAM Bounce Average Integration: Complete Implementation Guide

## Overview

This document provides a comprehensive explanation of the bounce average integration implementation in DREAM (Disruption Runaway Electron Analysis Model), combining theoretical foundations from the research paper with actual code implementation details.

---

## Table of Contents

1. [Theoretical Foundation](#1-theoretical-foundation)
2. [极向角积分的数值处理详解](#2-极向角积分的数值处理详解-poloidal-angle-integration)
3. [Code Architecture](#3-code-architecture)
4. [Implementation Details](#4-implementation-details)
5. [Numerical Integration Methods](#5-numerical-integration-methods)
6. [Special Cases and Singularities](#6-special-cases-and-singularities)
7. [Optimization Strategies](#7-optimization-strategies)
8. [Complete Integration Flow](#8-complete-integration-flow)
9. [Key Equations](#9-key-equations)
10. [Practical Considerations](#10-practical-considerations)
11. [Summary](#11-summary)

---

## 1. Theoretical Foundation

### 1.1 Bounce Average Definition

In DREAM, electrons are modeled using a **bounce-averaged drift-kinetic equation**:

$$\frac{\partial f}{\partial t} = \sum_{m,n} \frac{1}{\mathcal{V}'} \frac{\partial}{\partial z^m} \left[ \mathcal{V}' \left(-\{A^m\}f + \{D^{mn}\}\frac{\partial f}{\partial z^n}\right) \right] + \{S\}$$

where the **bounce average** is defined as:

$$\{X\} = \frac{1}{\mathcal{V}'} \int_0^{2\pi} d\zeta \int_0^{2\pi} d\phi \oint d\theta \sqrt{g} X$$

$$\mathcal{V}' = \int_0^{2\pi} d\zeta \int_0^{2\pi} d\phi \oint d\theta \sqrt{g}$$

### 1.2 Coordinate System

DREAM uses constants of motion as coordinates:
- **r**: Flux-surface label (distance from magnetic axis at B_min)
- **p**: Particle momentum magnitude (normalized to m_e c)
- **ξ₀**: Pitch angle at minimum magnetic field: ξ₀ = (B·p)/(Bp)|_{B=B_min}

### 1.3 Phase-Space Metric

The phase-space Jacobian is:

$$\sqrt{g} = p^2 \frac{B}{B_{\min}} \frac{\xi_0}{\xi} \mathcal{J}$$

where:
- $\mathcal{J} = \frac{1}{|\nabla\phi \cdot (\nabla\theta \times \nabla r)|}$ is the spatial Jacobian
- $\xi = \text{sgn}(\xi_0)\sqrt{1 - (1-\xi_0^2)B/B_{\min}}$ is the local pitch

### 1.4 Poloidal Orbit Integral

The poloidal integral depends on particle type:

$$\oint d\theta X = \begin{cases} 
\int_{-\pi}^{\pi} d\theta X, & |\xi_0| > \xi_T \quad \text{(passing)} \\
\int_{\theta_{b1}}^{\theta_{b2}} d\theta [X(\xi) + X(-\xi)], & 0 < \xi_0 \leq \xi_T \quad \text{(trapped)} \\
0, & -\xi_T \leq \xi_0 \leq 0 \quad \text{(negative trapped)}
\end{cases}$$

where:
- **Bounce points** θ_b1, θ_b2: angles where ξ(ξ₀, θ) = 0
- **Trapped-passing boundary**: ξ_T = √(1 - B_min/B_max)

### 1.5 The Three Integration Variables

The bounce average involves integration over **three angular variables**:

#### 1. Toroidal Angle (φ) → Analytically Integrated
- Integration range: ∫₀²π dφ
- Gives factor of **2π** (axisymmetric assumption)

#### 2. Gyrophase (ζ) → Analytically Integrated  
- Integration range: ∫₀²π dζ
- Gives factor of **2π** (zero orbit-width limit)

#### 3. Poloidal Angle (θ) → Requires Numerical Integration
- Integration range: ∮ dθ (depends on particle type)
- **This is the core computational challenge**
- Detailed numerical treatment in [Section 2](#2-极向角积分的数值处理详解-poloidal-angle-integration)

**Why only θ needs numerical integration?**

Under axisymmetric, zero orbit-width assumptions, the integrand becomes independent of φ and ζ:

$$\{X\} = \frac{1}{\mathcal{V}'} \underbrace{\int_0^{2\pi} d\zeta}_{2\pi} \underbrace{\int_0^{2\pi} d\phi}_{2\pi} \oint d\theta \sqrt{g} X = \frac{4\pi^2}{\mathcal{V}'} \oint d\theta \sqrt{g} X$$

In code, one 2π factor appears explicitly (BounceAverager.cpp line 598):
```cpp
BounceIntegral += 2*M_PI * w * Metric_data[it] * Function;
```
The other 2π cancels during normalization.

**The poloidal angle integration is complex because:**
1. Integration limits depend on particle type (passing vs trapped)
2. Bounce point finding required for trapped particles
3. Square-root singularities at bounce points
4. Magnetic geometry varies along the orbit

→ **See Section 2 for complete numerical treatment!**

---

## 2. 极向角积分的数值处理详解 (Poloidal Angle Integration)

本节详细说明DREAM中极向角θ的数值积分处理方法，这是bounce average实现的核心内容。

### 2.1 核心问题：积分区间依赖粒子类型

```
通行粒子 (passing):    θ ∈ [-π, π]         完整轨道
捕获粒子 (trapped):    θ ∈ [θ_b1, θ_b2]    两个反弹点之间
负pitch捕获:           积分为0              避免重复计算
```

### 2.2 第一步：判断粒子类型

```cpp
// BounceAverager.cpp line 421
bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir, i, j, fluxGridType, grid);

// 判断依据：比较 (1-ξ₀²) 和 B_min/B_max
real_t xi0 = GetXi0(ir, i, j, fluxGridType);
real_t Bmin = fluxSurfaceAverager->GetBmin(ir, fluxGridType);
real_t Bmax = fluxSurfaceAverager->GetBmax(ir, fluxGridType);

// 如果 (1-ξ₀²) > B_min/B_max，则是捕获粒子
if ((1 - xi0*xi0) > Bmin/Bmax) {
    isTrapped = true;
}
```

**物理意义**：
- ξ₀² = v_∥²/v² 在B_min处
- 如果粒子在B_max处的平行速度会变为虚数（即v_∥² < 0），说明粒子无法到达B_max，被"捕获"在磁阱中

### 2.3 第二步：寻找反弹点（仅捕获粒子需要）

对于捕获粒子，需要找到两个反弹点θ_b1和θ_b2，满足：

$$\xi(\theta_b) = \sqrt{1 - (1-\xi_0^2)\frac{B(\theta_b)}{B_{\min}}} = 0$$

等价于求解方程：

$$f(\theta) = 1 - (1-\xi_0^2)\frac{B(\theta)}{B_{\min}} = 0$$

**为什么需要数值找根？**

即使B(θ)有解析表达式（如Shafranov模型），上述方程仍然是**超越方程**，无法解析求解θ_b，必须使用数值方法。

**DREAM的实现方式**：

```cpp
// 调用 FindBouncePoints 函数
FluxSurfaceAverager::FindBouncePoints(
    ir,              // 径向位置
    Bmin,            // 最小磁场
    theta_Bmin,      // B_min 的θ位置
    theta_Bmax,      // B_max 的θ位置
    FSA,             // FluxSurfaceAverager指针
    xi0,             // pitch角
    fluxGridType,    // 网格类型
    &theta_b1,       // 输出：第一个反弹点
    &theta_b2,       // 输出：第二个反弹点
    gsl_fsolver,     // GSL求根器
    isSymmetric      // 是否对称磁场
);
```

**内部工作原理**（详见[4.2节](#42-bounce-point-finding)）：

1. **定义目标函数** f(θ) = ξ²(θ)
2. **确定搜索区间** [θ_Bmin, θ_Bmax]
3. **使用GSL Brent方法**分两步找根：
   - 第一步：在 [θ_Bmin, θ_Bmax] 找 θ_b2
   - 第二步：利用对称性或搜索另一半区间找 θ_b1
4. **返回两个反弹点** θ_b1, θ_b2

**关键依赖**：
- 需要知道 B(θ) 的值 → 通过 `GeometricQuantitiesAtTheta()` 获取
- B(θ) 的来源 → 解析模型或EFIT插值（见下文说明）
- 找根算法 → GSL Brent方法（详见[4.2节](#42-bounce-point-finding)）

**关键问题：B(θ) 的函数表达式是什么？**

在DREAM中，**B(θ) 没有简单的解析表达式**，而是通过以下方式获得：

1. **数值平衡数据**：
   - 磁场数据来自EFIT等平衡代码的输出
   - 以离散网格点形式存储：B(ir, θ_k)，k = 0, 1, ..., N_θ
   - 典型的角向分辨率：N_θ ≈ 100-200 个点

2. **插值方法**：
   ```cpp
   // AnalyticBRadialGridGenerator.cpp line 237
   void EvaluateGeometricQuantities(len_t ir, real_t theta, 
                                    real_t &B, real_t &Jacobian, ...) {
       // 对于解析模型（如CQL3D）：
       // B(θ) 通过Shafranov位移、椭圆度、三角形变等参数计算
       
       real_t ct = cos(theta);
       real_t st = sin(theta);
       
       // 计算 Shafranov 位移项
       real_t cdt = cos(delta[ir]*st);
       real_t sdt = sin(delta[ir]*st);
       
       // 计算 R/R0
       ROverR0 = 1 + (Delta[ir] + r[ir]*cos(theta+delta*sin(theta)))/R0;
       
       // 计算总磁场 B = sqrt(B_tor^2 + B_pol^2)
       real_t Btor = BtorGOverR0[ir]/ROverR0;  // 环向场
       real_t BpolSq = NablaR2 * (psiPrime/(2π*ROverR0))^2;  // 极向场
       B = sqrt(Btor*Btor + BpolSq);
   }
   ```

3. **两种实现方式**：
   
   **a) AnalyticBRadialGridGenerator**（解析模型 - 本例使用）：
   - 使用参数化平衡模型（Shafranov几何）
   - B(θ) 有半解析表达式，包含三角函数
   - 适用于简化物理研究
   
   **具体公式**（AnalyticBRadialGridGenerator.cpp line 237-272）：
   ```cpp
   void EvaluateGeometricQuantities(len_t ir, real_t theta, ...) {
       real_t ct = cos(theta);
       real_t st = sin(theta);
       
       // Shafranov 位移项
       real_t cdt = cos(delta[ir]*st);  // δ(r) = Shafranov位移
       real_t sdt = sin(delta[ir]*st);
       
       // 计算 R/R0 (大半径比)
       ROverR0 = 1 + (Delta[ir] + r[ir]*cos(theta+delta*sin(theta)))/R0;
       
       // 计算雅可比行列式
       Jacobian = r[ir] * ROverR0 * JOverRr;
       
       // 计算 |∇r|²
       NablaR2 = (kappaTerm² + deltaTerm²) / JOverRr²;
       
       // 总磁场 B = sqrt(B_tor² + B_pol²)
       real_t Btor = BtorGOverR0[ir]/ROverR0;          // 环向场 ~ 1/R
       real_t BpolSq = NablaR2 * (psiPrime/(2π*ROverR0))²;  // 极向场
       B = sqrt(Btor*Btor + BpolSq);
   }
   ```
   
   **关键参数**（来自Python设置）：
   ```python
   ds.radialgrid.setB0(B0)              # B0 = 1.4 T (轴心磁场)
   ds.radialgrid.setShaping(psi=0.0001, GOverR0=B0)  # G/R0 = B0 (环向场参数)
   ds.radialgrid.setMinorRadius(a)      # a = 0.01 m (小半径)
   ds.radialgrid.setMajorRadius(R)      # R0 = 1.67 m (大半径)
   ```
   
   这些参数定义了完整的Shafranov平衡，B(θ) 可以直接通过上述公式计算。
   
   **b) NumericBRadialGridGenerator**（数值平衡）：
   - 从EFIT g-file读取真实平衡数据
   - B(ir, θ) 存储在二维数组中
   - 通过**双线性插值**或**样条插值**获得任意θ的B值
   - 适用于真实托卡马克模拟

4. **找根过程中的B(θ)调用流程**：
   ```
   FindBouncePoints()
     └─> xiParticleFunction(theta)
          └─> GeometricQuantitiesAtTheta(ir, theta, B, ...)
               ├─> Analytic: 直接计算 B(θ) 公式
               └─> Numeric:  插值 B(ir, θ_k) 数据
   ```

5. **为什么需要数值找根？**
   
   即使对于解析模型，方程 `1 - (1-ξ₀²)B(θ)/B_min = 0` 也**无法解析求解**，因为：
   - B(θ) 包含复杂的三角函数组合
   - 涉及 Shafranov 位移 δ(r)、椭圆度 κ(r) 等径向依赖参数
   - 方程形式：`B(θ) = B_min/(1-ξ₀²)` 是超越方程
   
   因此必须使用**数值根查找算法**（如Brent方法）。

**找根策略**：

```cpp
// 第一个根在 [theta_Bmin, theta_Bmax] 区间
x_lower = theta_Bmin;
x_upper = theta_Bmax;
FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);
*theta2 = (gsl_func.function(x_lower) >= 0) ? x_lower : x_upper;

// 第二个根利用对称性或搜索另一半区间
if (isSymmetric) {
    *theta1 = -*theta2;  // 上下对称场，直接镜像
} else {
    // 非对称场，在另一半区间搜索
    x_upper = theta_Bmin;
    x_lower = theta_Bmax - 2*M_PI;
    FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);
    *theta1 = ...;
}
```

### 2.4 第三步：选择数值积分方法

**目标：计算极向角积分**

在找到反弹点 θ_b1, θ_b2 后，需要计算以下积分：

$$\text{BounceIntegral} = \int_{\theta_{b1}}^{\theta_{b2}} \sqrt{g}(\theta) \cdot X(\theta) \, d\theta$$

其中：
- **相空间度规**：$\sqrt{g}(\theta) = p^2 \frac{B(\theta)}{B_{\min}} \frac{\xi_0}{\xi(\theta)} \mathcal{J}(\theta)$
- **被积函数 X(θ)**：取决于具体物理量，例如：
  - 计算相空间体积：X = 1
  - 电场加速项：X = E_∥(θ) · ξ(θ)
  - 碰撞算符：X = ν_D(θ) · ξ²(θ)/ξ₀²
  - 同步辐射：X = B³(θ) 或 B²(θ)·ξ²(θ)/ξ₀²
  
  X(θ) 通常由幂次组合构成：X = (ξ/ξ₀)^a · (B/B_min)^b · (R/R₀)^c · |∇r|^d

**关键挑战**：
- 被积函数在反弹点有**平方根奇异性**：$\sqrt{g} \propto \frac{1}{\xi(\theta)} \propto \frac{1}{\sqrt{\theta_b - \theta}}$
- 需要选择合适的数值积分方法来处理奇异性

DREAM支持三种积分方法：

#### **方法A：固定节点Quadrature（默认）**

**Chebyshev quadrature**（捕获粒子默认）：

```cpp
// BounceAverager.cpp lines 116-120
case QUAD_FIXED_CHEBYSHEV:
    quadratureRule = gsl_integration_fixed_chebyshev;
    // 权重函数 w(x) = 1/√((x_max-x)(x-x_min))
    QuadFunc = [](real_t x, real_t xmin, real_t xmax) {
        return 1/sqrt((xmax-x)*(x-xmin));
    };
    break;
```

**为什么用Chebyshev？**
- 捕获粒子的相空间度规√g在反弹点有**平方根奇异性**：√g ∝ 1/√(θ_b - θ)
- Chebyshev权重函数恰好能抵消这个奇异性
- 比Legendre更准确，比adaptive更快

**实现过程**：

```cpp
// 1. 在参考区间[0,1]上生成quadrature点和权重
gsl_w = gsl_integration_fixed_alloc(quadratureRule, ntheta_interp_trapped, 
                                    0, 1, 0, 0);
theta_trapped_ref = gsl_w->x;       // 参考点 x_i ∈ [0,1]
weights_trapped_ref = gsl_w->weights; // 参考权重 w_i

// 2. 归一化权重（除以权重函数）
for(len_t it=0; it<ntheta_interp_trapped; it++)
    weights_trapped_ref[it] /= QuadFunc(theta_trapped_ref[it], 0, 1);

// 3. 映射到实际积分区间 [theta_b1, theta_b2]
weightScaleFactor = theta_b2 - theta_b1;
w = weightScaleFactor * weights_trapped_ref[it];

// 4. 在quadrature点上评估被积函数并求和
for (len_t it = 0; it < ntheta; it++) {
    // 将参考点映射到实际θ
    real_t theta = theta_b1 + theta_trapped_ref[it] * (theta_b2 - theta_b1);
    
    // 获取几何量在θ点的值（已预插值）
    real_t xiOverXi0 = MomentumGrid::evaluateXiOverXi0(xi0, BOverBmin_data[it]);
    real_t Function = AssembleBAFunc(xiOverXi0, BOverBmin_data[it],
                                    ROverR0_data[it], NablaR2_data[it],
                                    Flist_eval);
    
    // 累加：∑ w_i × √g(θ_i) × F(θ_i)
    BounceIntegral += 2*M_PI * w * Metric_data[it] * Function;
}
```

**Legendre quadrature**（通行粒子可用）：

```cpp
case QUAD_FIXED_LEGENDRE:
    quadratureRule = gsl_integration_fixed_legendre;
    QuadFunc = [](real_t x, real_t xmin, real_t xmax) { return 1; };
    break;
```

- 适用于光滑函数
- 权重函数w(x)=1，没有特殊处理奇异性

#### **方法B：自适应积分（Adaptive Quadrature）**

当被积函数复杂或精度要求高时使用：

```cpp
// BounceAverager.cpp lines 414-428
if ((!isTrapped && integratePassingAdaptive) || 
    (isTrapped && integrateTrappedAdaptive)) {
    
    // 检测是否有奇异性
    if (isTrapped && F_eval(0, 1, 1, 1, &params) != 0)
        params.integrateQAWS = true;  // 需要使用QAWS
    
    gsl_function GSL_func;
    GSL_func.function = &(FluxSurfaceAverager::BounceIntegralFunction);
    GSL_func.params = &params;
    
    real_t epsabs = 0, epsrel = 1e-6, lim = gsl_adaptive->limit, error;
    
    if (params.integrateQAWS)
        // QAWS：处理端点奇异性的自适应积分
        gsl_integration_qaws(&GSL_func, theta_b1, theta_b2, qaws_table,
                            epsabs, epsrel, lim, gsl_adaptive, 
                            &BounceIntegral, &error);
    else
        // QAG：一般自适应积分（Gauss-Kronrod规则）
        gsl_integration_qag(&GSL_func, theta_b1, theta_b2,
                           epsabs, epsrel, lim, QAG_KEY,
                           gsl_adaptive, &BounceIntegral, &error);
}
```

**QAG vs QAWS的区别**：

| 方法 | 适用场景 | 权重函数 |
|------|---------|---------|
| QAG | 光滑被积函数 | w(x) = 1 |
| QAWS | 端点有(x-a)^α(b-x)^β奇异性 | w(x) = (x-a)^α(b-x)^β |

**QAWS配置**：

```cpp
// BounceAverager.cpp line 65
qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
```

这对应权重函数：w(x) = (x-a)^{-1/2}(b-x)^{-1/2}，正好匹配捕获粒子的平方根奇异性。

**QAWS如何处理奇异性**：

```cpp
// FluxSurfaceAverager.cpp lines 339-358
real_t EvaluateBounceIntegrandWithQAWSWeightAtSingularity(...) {
    // 在反弹点附近做Taylor展开
    // ξ ≈ √[(dB/dθ)|_{θ_b} × (θ - θ_b) / B_min]
    
    real_t eps = 1e-4 * theta0;
    
    // 数值计算 dB/dθ 在反弹点
    real_t BPrimeAtThetaB = fabs(
        (FSA->BAtTheta(ir, theta0+eps, fgType) - 
         FSA->BAtTheta(ir, theta0-eps, fgType)) / (2*eps)
    );
    
    // Taylor展开因子：ξ/√(θ-θ_b) 的极限值
    real_t divergentSqrtFactor = 
        sqrt(xi0*xi0 / ((BPrimeAtThetaB/Bmin) * (1-xi0*xi0)));
    
    // 完整的被积函数 × √g
    real_t sqrtGTimesXiOverXi0 = 2*M_PI * B/Bmin;
    
    return sqrtGTimesXiOverXi0 * divergentSqrtFactor * otherSqrtFactor;
}
```

### 2.5 第四步：处理特殊情况

#### **情况1：捕获-通行边界附近的网格单元**

```cpp
// BounceAverager.cpp lines 434-436
if (fluxSurfaceAverager->shouldCellAverageBounceIntegral(ir, xi_f1, xi_f2, fluxGridType)) {
    // 网格单元跨越了捕获-通行边界，需要对ξ₀做cell average
    return fluxSurfaceAverager->EvaluateCellAveragedBounceIntegralOverP2(...);
}
```

**为什么要cell average？**
- 在ξ₀ = ξ_T处，相空间度规有对数奇异性
- 单个点的bounce average发散
- 需要对整个网格单元积分再平均

#### **情况2：ξ₀ = 0 的网格单元**

```cpp
// BounceAverager.cpp lines 439-445
else if (xi_f1 < 0 && xi_f2 > 0) {
    // 网格单元包含 ξ₀ = 0（最深捕获）
    fluxGridType = FLUXGRIDTYPE_P2;
    j += 1;
    isTrapped = BounceSurfaceQuantity::IsTrapped(ir, i, j, fluxGridType, grid);
    xi0 = GetXi0(ir, i, j, fluxGridType);
    SingularPointCorrection = xi0/2;  // 假设线性变化，取平均值
}
```

**物理意义**：
- ξ₀ = 0 时，θ_b1 = θ_b2 = θ_Bmin，反弹点重合
- bounce interval → 0，积分 → 0
- 假设被积函数从ξ₀=0到ξ₀线性变化，取平均值ξ₀/2

#### **情况3：负pitch捕获粒子**

```cpp
// BounceAverager.cpp lines 460-462
if (isTrapped) {
    if (xi0 < 100*realeps)  // ξ₀ < 0 的捕获粒子
        return 0;  // 不独立存在，由正pitch描述
}
```

**原因**：
- 负pitch捕获粒子和正pitch是同一个物理轨道的两个方向
- 在正pitch的积分中已经通过 `[X(ξ) + X(-ξ)]` 包含了两个方向
- 避免重复计算

### 2.6 第五步：对称性优化

对于上下对称的托卡马克磁场：

```cpp
// FluxSurfaceAverager.cpp lines 202-205
if (geometryIsSymmetric)
    theta_max = M_PI;      // 只积分 0 到 π
else
    theta_max = 2*M_PI;    // 积分 0 到 2π

// 权重乘以2补偿
if (geometryIsSymmetric)
    for (len_t it=0; it<ntheta_interp; it++)
        weights[it] *= 2;
```

**捕获粒子的对称性**：

```cpp
// BounceAverager.cpp line 483-484
if (geometryIsSymmetric && theta_b1 > 0)
    SingularPointCorrection *= 2;  // 两个对称区间贡献相同
```

### 2.7 完整积分公式

最终数值积分计算的是：

$$\text{BounceIntegral} = 2\pi \sum_{i=1}^{N_\theta} w_i \cdot \sqrt{g}(\theta_i) \cdot X(\theta_i)$$

其中：
- $w_i$：quadrature权重（已包含区间缩放因子）
- $\sqrt{g}(\theta_i) = p^2 \frac{B(\theta_i)}{B_{\min}} \frac{\xi_0}{\xi(\theta_i)} \mathcal{J}(\theta_i)$
- $X(\theta_i)$：被积函数在θ_i的值
- $2\pi$：来自环向或回旋角的解析积分

**归一化得到bounce average**：

$$\{X\} = \frac{\text{BounceIntegral}}{\mathcal{V}'} = \frac{\int \sqrt{g} X d\theta}{\int \sqrt{g} d\theta}$$

### 2.8 总结：数值处理的关键点

1. **粒子分类** → 决定积分区间
2. **找反弹点** → Brent根查找算法
3. **选quadrature** → Chebyshev处理奇异性，Legendre用于光滑函数，Adaptive用于高精度
4. **处理奇异性** → QAWS或cell averaging
5. **利用对称性** → 减少一半计算量
6. **预插值几何量** → 提高计算效率

这套方法既保证了数值精度，又通过优化策略保持了计算效率！

---

## 3. Code Architecture

The bounce average implementation consists of three main components:

### 3.1 Grid.cpp (`/data/zhzhou/DREAM/fvm/Grid/Grid.cpp`)

**Role**: Main orchestration and interface

**Key responsibilities**:
- Creates and manages `BounceAverager` instance
- Calls `RebuildBounceAveragedQuantities()` to pre-compute bounce averages
- Provides public interface methods

**Key methods**:
```cpp
// Constructor - creates BounceAverager
Grid::Grid(RadialGrid *rg, MomentumGrid *mg, ...) {
    bounceAverager = new BounceAverager(this, FSA, ntheta_interp_trapped, qm_trapped);
}

// Rebuild all bounce-averaged quantities
void Grid::RebuildBounceAveragedQuantities() {
    SetBounceAveragePXi(BA_xi_fr, FLUXGRIDTYPE_RADIAL, ...);
    SetBounceAveragePXi(BA_xi_f1, FLUXGRIDTYPE_P1, ...);
    // ... compute various bounce averages
}

// Public interface for calculating bounce average
real_t Grid::CalculateBounceAverage(len_t ir, len_t i, len_t j, ...) {
    return bounceAverager->CalculateBounceAverage(ir, i, j, ...);
}
```

### 3.2 BounceAverager.cpp (`/data/zhzhou/DREAM/fvm/Grid/BounceAverager.cpp`)

**Role**: Core bounce averaging logic

**Key responsibilities**:
- Handles momentum-grid dependent calculations
- Manages quadrature rules for trapped/passing particles
- Computes bounce integrals using GSL numerical integration
- Determines trapped vs passing particle status

**Key data structures**:
```cpp
class BounceAverager {
private:
    Grid *grid;
    FluxSurfaceAverager *fluxSurfaceAverager;
    
    // Quadrature parameters
    len_t ntheta_interp_trapped;      // Number of points for trapped orbits
    len_t ntheta_interp_passing;      // Number of points for passing orbits
    real_t *theta_trapped_ref;        // Reference quadrature points
    real_t *weights_trapped_ref;      // Reference quadrature weights
    
    // Integration flags
    bool integrateTrappedAdaptive;    // Use adaptive quad for trapped?
    bool integratePassingAdaptive;    // Use adaptive quad for passing?
    
    // GSL integration workspaces
    gsl_integration_fixed_workspace *gsl_w;
    gsl_integration_workspace *gsl_adaptive;
    gsl_root_fsolver *gsl_fsolver;
    gsl_integration_qaws_table *qaws_table;
    
    // Pre-computed geometric quantities
    BounceSurfaceQuantity *BOverBmin;
    BounceSurfaceQuantity *ROverR0;
    BounceSurfaceQuantity *NablaR2;
    BounceSurfaceMetric *Metric;
    
    // Trapped particle data
    bool **isTrapped;                 // Is particle trapped?
    real_t **theta_b1, **theta_b2;   // Bounce points
};
```

**Key methods**:
```cpp
// Main calculation method
real_t CalculateBounceAverage(len_t ir, len_t i, len_t j, 
                              fluxGridType, function F, ...);

// Core bounce integral evaluation
real_t EvaluateBounceIntegralOverP2(len_t ir, len_t i, len_t j, ...);

// Initialize quadrature rules
void InitializeQuadrature(quadrature_method q_method);

// Find bounce points and trapped status
bool InitializeBounceIntegralQuantities();
```

### 3.3 FluxSurfaceAverager.cpp (`/data/zhzhou/DREAM/fvm/Grid/FluxSurfaceAverager.cpp`)

**Role**: Geometric quantities and flux surface operations

**Key responsibilities**:
- Magnetic field geometry interpolation
- Finding bounce points using root-finding
- Flux surface averaging
- Handling singular integrands

**Key methods**:
```cpp
// Find bounce points where ξ = 0
void FindBouncePoints(len_t ir, real_t Bmin, real_t theta_Bmin, 
                      real_t theta_Bmax, ..., real_t xi0, ...);

// Evaluate bounce integral at arbitrary (p, ξ₀)
real_t EvaluatePXiBounceIntegralAtP(len_t ir, real_t xi0, ...);

// Assemble bounce average function from powers
real_t AssembleBAFunc(real_t xiOverXi0, real_t BOverBmin, 
                      real_t ROverR0, real_t NablaR2, const int_t *Flist);

// Handle singular integrand at bounce points
real_t BounceIntegralFunction(real_t theta, void *par);
```

---

## 4. Implementation Details

### 4.1 Quadrature Methods

DREAM supports multiple integration methods for different scenarios:

#### Fixed Quadrature Rules

**Legendre Quadrature** (`QUAD_FIXED_LEGENDRE`):
- Best for smooth functions on finite intervals
- Weight function: w(x) = 1
- Used for passing particles or smooth integrands

**Chebyshev Quadrature** (`QUAD_FIXED_CHEBYSHEV`):
- Optimal for trapped orbits with square-root singularities
- Weight function: w(x) = 1/√((x_max-x)(x-x_min))
- Default choice for trapped particles

```cpp
// From BounceAverager.cpp lines 101-141
void BounceAverager::InitializeQuadrature(quadrature_method q_method) {
    const gsl_integration_fixed_type *quadratureRule;
    function<real_t(real_t,real_t,real_t)> QuadFunc;
    
    switch(q_method) {
        case QUAD_FIXED_LEGENDRE:
            quadratureRule = gsl_integration_fixed_legendre;
            QuadFunc = [](real_t x, real_t xmin, real_t xmax) { return 1; };
            break;
            
        case QUAD_FIXED_CHEBYSHEV:
            quadratureRule = gsl_integration_fixed_chebyshev;
            QuadFunc = [](real_t x, real_t xmin, real_t xmax) {
                return 1/sqrt((xmax-x)*(x-xmin));
            };
            break;
            
        case QUAD_ADAPTIVE:
            integrateTrappedAdaptive = true;
            return;
    }
    
    // Generate reference quadrature on [0,1]
    gsl_w = gsl_integration_fixed_alloc(quadratureRule, ntheta_interp_trapped, 
                                        0, 1, 0, 0);
    theta_trapped_ref = gsl_w->x;
    weights_trapped_ref = gsl_w->weights;
    
    // Normalize weights by dividing by quadrature weight function
    for(len_t it=0; it<ntheta_interp_trapped; it++)
        weights_trapped_ref[it] /= QuadFunc(theta_trapped_ref[it], 0, 1);
}
```

#### Adaptive Quadrature

**GSL QAG** (General adaptive quadrature):
- For regular integrands
- Uses Gauss-Kronrod rules

**GSL QAWS** (Adaptive quadrature for singular functions):
- Specifically for integrable singularities at endpoints
- Essential for trapped particles near bounce points
- Weight function: w(x) = (x-a)^α (b-x)^β ln^μ(x-a) ln^ν(b-x)

```cpp
// From BounceAverager.cpp lines 414-428
if (( !isTrapped && integratePassingAdaptive) || 
    (isTrapped && integrateTrappedAdaptive)) {
    
    if (isTrapped && F_eval(0,1,1,1,&params) != 0)
        params.integrateQAWS = true;  // Detect singularity
    
    gsl_function GSL_func;
    GSL_func.function = &(FluxSurfaceAverager::BounceIntegralFunction);
    GSL_func.params = &params;
    
    real_t epsabs = 0, epsrel = 1e-6, lim = gsl_adaptive->limit, error;
    
    if (params.integrateQAWS)
        gsl_integration_qaws(&GSL_func, theta_b1, theta_b2, qaws_table,
                            epsabs, epsrel, lim, gsl_adaptive, 
                            &BounceIntegral, &error);
    else
        gsl_integration_qag(&GSL_func, theta_b1, theta_b2,
                           epsabs, epsrel, lim, QAG_KEY,
                           gsl_adaptive, &BounceIntegral, &error);
}
```

### 4.2 Bounce Point Finding

Bounce points are found by solving for where the parallel velocity vanishes:

$$\xi(\xi_0, \theta) = \sqrt{1 - (1-\xi_0^2)\frac{B(\theta)}{B_{\min}}} = 0$$

This is equivalent to finding roots of:

$$\xi^2 = 1 - (1-\xi_0^2)\frac{B(\theta)}{B_{\min}} = 0$$

**Implementation using GSL Brent's method**:

```cpp
// From FluxSurfaceAverager.cpp lines 596-608
void FluxSurfaceAverager::FindBouncePoints(
    len_t ir, real_t Bmin, real_t theta_Bmin, real_t theta_Bmax,
    FluxSurfaceAverager *FSA, real_t xi0, fluxGridType fluxGridType,
    real_t *theta_b1, real_t *theta_b2, 
    gsl_root_fsolver* gsl_fsolver, bool isSymmetric) {
    
    // Define function ξ²(θ) to find roots
    xiFuncParams xi_params = {xi0, ir, Bmin, FSA, fluxGridType};
    gsl_function gsl_func;
    gsl_func.function = &(xiParticleFunction);  // Returns ξ²
    gsl_func.params = &xi_params;
    
    // Find two roots bracketing the magnetic well
    FindThetas(theta_Bmin, theta_Bmax, theta_b1, theta_b2,
               gsl_func, gsl_fsolver, isSymmetric);
}

// The function whose roots we seek
real_t FluxSurfaceAverager::xiParticleFunction(real_t theta, void *p) {
    struct xiFuncParams *params = (struct xiFuncParams *) p;
    real_t B, Jacobian, ROverR0, NablaR2;
    
    params->FSA->GeometricQuantitiesAtTheta(params->ir, theta, 
                                            B, Jacobian, ROverR0, NablaR2, 
                                            params->fgType);
    
    real_t Bmin = params->Bmin;
    real_t BOverBmin = (Bmin != 0) ? B/Bmin : 1;
    real_t xi0 = params->xi0;
    
    // Return ξ² = 1 - (1-ξ₀²)B/B_min
    return 1 - (1 - xi0*xi0) * BOverBmin;
}
```

**Root-finding strategy**:

```cpp
// From FluxSurfaceAverager.cpp lines 619-672
void FluxSurfaceAverager::FindThetas(
    real_t theta_Bmin, real_t theta_Bmax, 
    real_t *theta1, real_t *theta2,
    gsl_function gsl_func, gsl_root_fsolver *gsl_fsolver, 
    bool isSymmetric) {
    
    real_t root = 0;
    real_t x_lower = theta_Bmin, x_upper = theta_Bmax;
    
    // Adjust interval if Bmin is in lower half plane
    if (x_lower > x_upper)
        x_lower -= 2*M_PI;
    
    // Find first root (upper bounce point)
    FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);
    
    if (gsl_func.function(x_lower, gsl_func.params) >= 0)
        *theta2 = x_lower;
    else if (gsl_func.function(x_upper, gsl_func.params) >= 0)
        *theta2 = x_upper;
    else
        throw FVMException("Unable to find valid theta root.");
    
    if (isSymmetric) {
        // For symmetric fields, use mirror symmetry
        if (theta_Bmin > 0 && gsl_func.function(0, gsl_func.params) < 0) {
            // Special case: look in [0, theta_Bmin]
            x_lower = 0;
            x_upper = theta_Bmin;
        } else {
            *theta1 = -*theta2;  // Mirror symmetry
            return;
        }
    } else {
        // Asymmetric field: search remaining interval
        x_upper = theta_Bmin;
        x_lower = theta_Bmax - (theta_Bmin < theta_Bmax ? 2*M_PI : 0);
    }
    
    // Find second root (lower bounce point)
    FindRoot(&x_lower, &x_upper, &root, gsl_func, gsl_fsolver);
    
    if (gsl_func.function(x_lower, gsl_func.params) >= 0)
        *theta1 = x_lower;
    else if (gsl_func.function(x_upper, gsl_func.params) >= 0)
        *theta1 = x_upper;
    else
        throw FVMException("Unable to find valid theta root.");
}
```

### 4.3 Bounce Integral Evaluation

The core function evaluates:

$$\text{BounceIntegral}(X) = \int \sqrt{g} X \, d\phi \, d\theta \, d\zeta$$

```cpp
// From BounceAverager.cpp lines 316-461
real_t BounceAverager::EvaluateBounceIntegralOverP2(
    len_t ir, len_t i, len_t j, fluxGridType fluxGridType,
    real_t(*F_ref)(real_t,real_t,real_t,real_t,void*), 
    void *F_ref_par, const int_t *Flist) {
    
    // Special case: single pitch grid point
    if (np2[0] == 1) {
        if (Flist != nullptr) {
            if (Flist[0]%2 == 1)  // Odd function → zero
                return 0;
            int_t Flist_new[4] = {Flist[1], Flist[2], Flist[3], Flist[4]};
            return 2*M_PI/(Flist[0]+1) * 
                   fluxSurfaceAverager->EvaluateFluxSurfaceIntegral(...);
        }
    }
    
    // Check if cell should be mirrored (negative-pitch trapped)
    if (fluxGridType != FLUXGRIDTYPE_P2) {
        real_t xi_f2 = GetXi0(ir, i, j+1, FLUXGRIDTYPE_P2);
        real_t xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
        
        if (xi_f2 < 100*realeps && xi_f2 > -xiT && xiT > 0)
            return 0;  // Negative-pitch trapped → skip
    }
    
    // Get particle parameters
    real_t xi0 = GetXi0(ir, i, j, fluxGridType);
    bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir, i, j, fluxGridType, grid);
    real_t Bmin = fluxSurfaceAverager->GetBmin(ir, fluxGridType);
    real_t Bmax = fluxSurfaceAverager->GetBmax(ir, fluxGridType);
    
    real_t SingularPointCorrection = 1;
    
    // Handle special cases near singularities
    if (Bmin != Bmax && fluxGridType == FLUXGRIDTYPE_DISTRIBUTION) {
        real_t xi_f1 = GetXi0(ir, i, j, FLUXGRIDTYPE_P2);
        real_t xi_f2 = GetXi0(ir, i, j+1, FLUXGRIDTYPE_P2);
        
        if (xi_f1 > xi_f2) swap(xi_f1, xi_f2);
        
        // Cell contains trapped-passing boundary
        if (fluxSurfaceAverager->shouldCellAverageBounceIntegral(ir, xi_f1, xi_f2, fluxGridType)) {
            return fluxSurfaceAverager->EvaluateCellAveragedBounceIntegralOverP2(...);
        }
        // Cell contains ξ₀ = 0
        else if (xi_f1 < 0 && xi_f2 > 0) {
            fluxGridType = FLUXGRIDTYPE_P2;
            j += 1;
            isTrapped = BounceSurfaceQuantity::IsTrapped(ir, i, j, fluxGridType, grid);
            xi0 = GetXi0(ir, i, j, fluxGridType);
            SingularPointCorrection = xi0/2;  // Linear interpolation assumption
        }
    }
    
    // Prepare function list for power-law assembly
    int_t *Flist_eval = nullptr;
    int_t Flist_copy[5];
    if (Flist != nullptr) {
        for (len_t k=0; k<5; k++)
            Flist_copy[k] = Flist[k];
        Flist_eval = Flist_copy;
    }
    
    real_t(*F_eval)(real_t,real_t,real_t,real_t,void*);
    
    if (isTrapped) {
        // Negative-pitch trapped particles don't exist independently
        if (xi0 < 100*realeps)
            return 0;
        
        // Handle symmetry: odd functions vanish, even functions double
        if (Flist != nullptr) {
            if (Flist[0]%2 == 1)
                return 0;  // Odd in ξ → zero
            else
                Flist_eval[4] *= 2;  // Even in ξ → multiply by 2
        }
        
        // Sum over both directions for trapped particles
        F_eval = FluxSurfaceAverager::BA_FUNC_TRAPPED;
    } else {
        F_eval = FluxSurfaceAverager::BA_FUNC_PASSING;
    }
    
    // Get bounce points
    real_t theta_b1 = BounceSurfaceQuantity::Theta_B1(ir, i, j, fluxGridType, grid);
    real_t theta_b2 = BounceSurfaceQuantity::Theta_B2(ir, i, j, fluxGridType, grid);
    
    // Account for symmetric geometry (two identical intervals)
    if (geometryIsSymmetric && theta_b1 > 0)
        SingularPointCorrection *= 2;
    
    // Set up integration parameters
    real_t BounceIntegral = 0;
    FluxSurfaceAverager::BounceIntegralParams params = {
        ir, xi0, theta_b1, theta_b2, fluxGridType, Bmin,
        F_ref, F_eval, F_ref_par, Flist_eval, fluxSurfaceAverager, false
    };
    
    // Choose integration method
    if ((!isTrapped && integratePassingAdaptive) || 
        (isTrapped && integrateTrappedAdaptive)) {
        
        // Check if integrand is singular (for QAWS)
        if (isTrapped && F_eval(0, 1, 1, 1, &params) != 0)
            params.integrateQAWS = true;
        
        gsl_function GSL_func;
        GSL_func.function = &(FluxSurfaceAverager::BounceIntegralFunction);
        GSL_func.params = &params;
        
        real_t epsabs = 0, epsrel = 1e-6, lim = gsl_adaptive->limit, error;
        
        if (params.integrateQAWS)
            gsl_integration_qaws(&GSL_func, theta_b1, theta_b2, qaws_table,
                                epsabs, epsrel, lim, gsl_adaptive,
                                &BounceIntegral, &error);
        else
            gsl_integration_qag(&GSL_func, theta_b1, theta_b2,
                               epsabs, epsrel, lim, QAG_KEY,
                               gsl_adaptive, &BounceIntegral, &error);
        
        return SingularPointCorrection * BounceIntegral;
    }
    
    // Fixed quadrature: use pre-computed values
    const real_t *BOverBmin_data = this->BOverBmin->GetData(ir, i, j, fluxGridType);
    const real_t *ROverR0_data = this->ROverR0->GetData(ir, i, j, fluxGridType);
    const real_t *NablaR2_data = this->NablaR2->GetData(ir, i, j, fluxGridType);
    const real_t *Metric_data = this->Metric->GetData(ir, i, j, fluxGridType);
    
    len_t ntheta;
    const real_t *weights;
    real_t weightScaleFactor;
    
    if (isTrapped) {
        ntheta = ntheta_interp_trapped;
        weights = weights_trapped_ref;
        weightScaleFactor = theta_b2 - theta_b1;
    } else {
        ntheta = ntheta_interp_passing;
        weights = this->weights_passing;
        weightScaleFactor = 1;
    }
    
    bool hasFList = (Flist_eval != nullptr);
    
    for (len_t it = 0; it < ntheta; it++) {
        real_t xiOverXi0 = MomentumGrid::evaluateXiOverXi0(xi0, BOverBmin_data[it]);
        real_t w = weightScaleFactor * weights[it];
        
        real_t Function = hasFList ?
            fluxSurfaceAverager->AssembleBAFunc(xiOverXi0, BOverBmin_data[it],
                                               ROverR0_data[it], NablaR2_data[it],
                                               Flist_eval)
            : F_eval(xiOverXi0, BOverBmin_data[it], ROverR0_data[it],
                    NablaR2_data[it], &params);
        
        BounceIntegral += 2*M_PI * w * Metric_data[it] * Function;
    }
    
    return SingularPointCorrection * BounceIntegral;
}
```

### 4.4 Calculating the Bounce Average

The final bounce average is computed by normalizing the integral:

```cpp
// From BounceAverager.cpp lines 280-306
real_t BounceAverager::CalculateBounceAverage(
    len_t ir, len_t i, len_t j, fluxGridType fluxGridType,
    real_t(*F)(real_t,real_t,real_t,real_t,void*), 
    void *par, const int_t *F_list) {
    
    real_t Vp, p, preFactor;
    
    // Get momentum at this grid point
    if (fluxGridType == FLUXGRIDTYPE_P1)
        p = grid->GetMomentumGrid(0)->GetP_f1(i, j);
    else if (fluxGridType == FLUXGRIDTYPE_P2)
        p = grid->GetMomentumGrid(0)->GetP_f2(i, j);
    else
        p = grid->GetMomentumGrid(0)->GetP(i, j);
    
    // Handle p=0 singularity
    if (p == 0) {
        Vp = grid->GetVpOverP2AtZero(ir)[j];
        preFactor = 1.0;
    } else {
        Vp = GetVp(ir, i, j, fluxGridType);
        preFactor = p*p;
    }
    
    if (!Vp)
        return 0;
    
    // Evaluate bounce integral
    real_t BI = preFactor * EvaluateBounceIntegralOverP2(
        ir, i, j, fluxGridType, F, par, F_list);
    
    if (!BI)
        return 0;
    
    // Normalize by phase-space volume
    return BI / Vp;
}
```

---

## 5. Numerical Integration Methods

### 5.1 Integration Strategy Selection

The choice of integration method depends on:

1. **Particle type**: trapped vs passing
2. **Integrand smoothness**: regular vs singular
3. **Performance requirements**: speed vs accuracy
4. **Grid resolution**: fine vs coarse

### 5.2 Fixed Quadrature Implementation

For fixed quadrature, the integral is approximated as:

$$\int_a^b w(x) f(x) dx \approx \sum_{i=1}^N w_i f(x_i)$$

where the reference quadrature on [0,1] is mapped to [θ_b1, θ_b2]:

```cpp
// Mapping from reference interval [0,1] to [theta_b1, theta_b2]
weightScaleFactor = theta_b2 - theta_b1;
w = weightScaleFactor * weights_ref[it];
```

### 5.3 Adaptive Quadrature Implementation

Adaptive quadrature automatically refines the integration grid based on error estimates:

```cpp
// Parameters for adaptive integration
real_t epsabs = 0;          // Absolute error tolerance
real_t epsrel = 1e-6;       // Relative error tolerance
len_t lim = gsl_adaptive->limit;  // Maximum subdivisions
real_t error;               // Estimated error

// General adaptive quadrature (QAG)
gsl_integration_qag(&GSL_func, theta_b1, theta_b2,
                   epsabs, epsrel, lim, 
                   GSL_INTEG_GAUSS41,  // 41-point Gauss-Kronrod rule
                   gsl_adaptive, 
                   &result, &error);

// Adaptive quadrature for singular integrands (QAWS)
gsl_integration_qaws(&GSL_func, theta_b1, theta_b2, 
                    qaws_table,  // Singularity parameters
                    epsabs, epsrel, lim,
                    gsl_adaptive,
                    &result, &error);
```

### 5.4 QAWS Singularity Handling

For trapped particles, the integrand has square-root singularities at bounce points. The QAWS table is configured with:

```cpp
// From BounceAverager.cpp line 65
qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
```

This corresponds to weight function: w(x) = (x-a)^{-1/2} (b-x)^{-1/2}

The integrand is rescaled to remove the singularity:

```cpp
// From FluxSurfaceAverager.cpp lines 339-358
real_t EvaluateBounceIntegrandWithQAWSWeightAtSingularity(
    len_t ir, real_t xi0, real_t theta, 
    real_t theta_b1, real_t theta_b2,
    real_t B, real_t Bmin, 
    FluxSurfaceAverager *FSA, fluxGridType fgType) {
    
    real_t otherSqrtFactor, theta0;
    
    // Determine which bounce point we're near
    if (theta_b2 - theta < theta - theta_b1) {
        otherSqrtFactor = sqrt(theta - theta_b1);
        theta0 = theta_b2;
    } else {
        otherSqrtFactor = sqrt(theta_b2 - theta);
        theta0 = theta_b1;
    }
    
    real_t eps = 1e-4 * theta0;
    
    // Compute dB/dθ at bounce point using finite difference
    real_t BPrimeAtThetaB = fabs(
        (FSA->BAtTheta(ir, theta0+eps, fgType) - 
         FSA->BAtTheta(ir, theta0-eps, fgType)) / (2*eps)
    );
    
    // Taylor expansion factor
    real_t divergentSqrtFactorOverXiOverXi0 = 
        sqrt(xi0*xi0 / ((BPrimeAtThetaB/Bmin) * (1-xi0*xi0)));
    
    // Explicit expression for metric × (ξ/ξ₀)
    real_t sqrtGTimesXiOverXi0 = 2*M_PI * B/Bmin;
    
    return sqrtGTimesXiOverXi0 * divergentSqrtFactorOverXiOverXi0 * otherSqrtFactor;
}
```

---

## 6. Special Cases and Singularities

### 6.1 Trapped-Passing Boundary

At the trapped-passing boundary ξ₀ = ξ_T, the phase-space metric has a logarithmic singularity. DREAM handles this through **cell averaging**:

```cpp
// Check if cell contains the trapped-passing boundary
if (fluxSurfaceAverager->shouldCellAverageBounceIntegral(ir, xi_f1, xi_f2, fluxGridType)) {
    return fluxSurfaceAverager->EvaluateCellAveragedBounceIntegralOverP2(
        ir, xi_f1, xi_f2, fluxGridType, F_ref, F_ref_par, Flist);
}
```

The cell average integrates over the pitch range within the cell:

$$\langle X \rangle_{\text{cell}} = \frac{1}{\Delta\xi_0} \int_{\xi_{0,1}}^{\xi_{0,2}} \{X\}(\xi_0) d\xi_0$$

### 6.2 Zero Momentum (p = 0)

At p = 0, the metric √g ∝ p² vanishes. Special handling:

```cpp
if (p == 0) {
    Vp = grid->GetVpOverP2AtZero(ir)[j];  // Pre-computed Vp/p² at p=0
    preFactor = 1.0;
} else {
    Vp = GetVp(ir, i, j, fluxGridType);
    preFactor = p*p;
}
```

### 6.3 Zero Pitch (ξ₀ = 0)

At ξ₀ = 0 (infinitely deeply trapped particles), the bounce points coincide: θ_b1 = θ_b2 = θ_Bmin

```cpp
if (std::abs(xi0) < 100*realeps) {  // ξ₀ ≈ 0
    theta_b1[ir][pind] = theta_Bmin;
    theta_b2[ir][pind] = theta_Bmin;
    // Integral goes to zero as bounce interval vanishes
}
```

For cells containing ξ₀ = 0, linear interpolation is assumed:

```cpp
else if (xi_f1 < 0 && xi_f2 > 0) {
    // Cell spans ξ₀ = 0
    SingularPointCorrection = xi0/2;  // Average over cell assuming linear variation
}
```

### 6.4 Negative-Pitch Trapped Particles

Negative-pitch trapped particles (-ξ_T ≤ ξ₀ ≤ 0) are not independent; their dynamics are described by positive-pitch counterparts. To avoid double counting:

```cpp
// From BounceAverager.cpp lines 329-334
if (fluxGridType != FLUXGRIDTYPE_P2) {
    real_t xi_f2 = GetXi0(ir, i, j+1, FLUXGRIDTYPE_P2);
    real_t xiT = grid->GetRadialGrid()->GetXi0TrappedBoundary(ir);
    
    if (xi_f2 < 100*realeps && xi_f2 > -xiT && xiT > 0)
        return 0;  // Skip negative-pitch trapped cells
}
```

### 6.5 Symmetric Geometry

For tokamaks with up-down symmetric magnetic fields, only half the orbit needs to be integrated:

```cpp
// From FluxSurfaceAverager.cpp lines 202-205
if (geometryIsSymmetric)
    theta_max = M_PI;      // Integrate from 0 to π
else
    theta_max = 2*M_PI;    // Integrate from 0 to 2π

// Multiply result by 2 for symmetric case
if (geometryIsSymmetric)
    for (len_t it=0; it<ntheta_interp; it++)
        weights[it] *= 2;
```

For trapped particles in symmetric fields:

```cpp
if (geometryIsSymmetric && theta_b1 > 0)
    SingularPointCorrection *= 2;  // Two identical intervals contribute
```

---

## 7. Optimization Strategies

### 7.1 PXi Grid Optimization

For p-ξ grids, the phase-space volume scales as Vp ∝ p². This allows significant optimization:

```cpp
// From BounceAverager.cpp lines 250-272
void BounceAverager::SetVpsPXi(real_t**&Vp, real_t**&Vp_f1, real_t**&VpOverP2) {
    Vp = new real_t*[nr];
    Vp_f1 = new real_t*[nr];
    VpOverP2 = new real_t*[nr];
    
    const real_t *p = grid->GetMomentumGrid(0)->GetP1();
    const real_t *p_f = grid->GetMomentumGrid(0)->GetP1_f();
    
    for (len_t ir = 0; ir < nr; ir++) {
        len_t n1 = np1[ir];
        len_t n2 = np2[ir];
        
        Vp[ir] = new real_t[n1*n2];
        Vp_f1[ir] = new real_t[(n1+1)*n2];
        VpOverP2[ir] = new real_t[n2];
        
        for (len_t j = 0; j < n2; j++) {
            // Compute once per pitch angle (independent of p)
            VpOverP2[ir][j] = EvaluateBounceIntegralOverP2(
                ir, 0, j, FLUXGRIDTYPE_DISTRIBUTION,
                RadialGrid::BA_FUNC_UNITY, nullptr, 
                RadialGrid::BA_PARAM_UNITY
            );
            
            // Scale by p² for all momentum points
            for (len_t i = 0; i < n1; i++) {
                Vp[ir][j*n1+i] = VpOverP2[ir][j] * p[i]*p[i];
                Vp_f1[ir][j*(n1+1)+i] = VpOverP2[ir][j] * p_f[i]*p_f[i];
            }
            Vp_f1[ir][j*(n1+1)+n1] = VpOverP2[ir][j] * p_f[n1]*p_f[n1];
        }
    }
}
```

**Performance benefit**: Reduces computation from O(N_p × N_ξ) to O(N_ξ) for Vp calculation.

### 7.2 Pre-computation of Geometric Quantities

Common geometric factors are pre-computed and stored:

```cpp
// From BounceAverager.cpp lines 42-45
BOverBmin = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetBOverBmin());
ROverR0 = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetROverR0());
NablaR2 = new BounceSurfaceQuantity(grid, fluxSurfaceAverager->GetNablaR2());
Metric = new BounceSurfaceMetric(grid, fluxSurfaceAverager->GetJacobian(),
                                 fluxSurfaceAverager->GetBOverBmin(),
                                 fluxSurfaceAverager);
```

These are interpolated to quadrature points during `Rebuild()`:

```cpp
// From BounceAverager.cpp lines 183-188
if (storeTrapped) {
    BOverBmin->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
    ROverR0->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
    NablaR2->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
    Metric->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
}
```

### 7.3 Power-Law Function Assembly

For common functions that are products of powers, DREAM uses an efficient assembly method:

```cpp
// From FluxSurfaceAverager.cpp lines 514-541
real_t FluxSurfaceAverager::AssembleBAFunc(
    real_t xiOverXi0, real_t BOverBmin, 
    real_t ROverR0, real_t NablaR2, 
    const int_t *Flist) {
    
    // Flist = [power_xi, power_B, power_R, power_nablaR2, coefficient]
    real_t BA_Func = Flist[4];  // Start with coefficient
    
    // Multiply by powers (positive and negative)
    if (Flist[0] > 0)
        for (int_t k=0; k<Flist[0]; k++)
            BA_Func *= xiOverXi0;
    else if (Flist[0] < 0)
        for (int_t k=0; k<-Flist[0]; k++)
            BA_Func /= xiOverXi0;
    
    if (Flist[1] > 0)
        for (int_t k=0; k<Flist[1]; k++)
            BA_Func *= BOverBmin;
    else if (Flist[1] < 0)
        for (int_t k=0; k<-Flist[1]; k++)
            BA_Func /= BOverBmin;
    
    // ... similar for ROverR0 and NablaR2
    
    return BA_Func;
}
```

This avoids expensive function evaluations for terms like:
- ξ/ξ₀, (ξ/ξ₀)², etc.
- B/B_min, (B/B_min)², etc.
- Mixed products like (ξ/ξ₀)²(B/B_min)³

### 7.4 Cached Bounce Averages

Commonly used bounce averages are pre-computed and stored:

```cpp
// From Grid.cpp lines 259-267
SetBounceAveragePXi(BA_xi_fr, FLUXGRIDTYPE_RADIAL, 
                   RadialGrid::BA_FUNC_XI, nullptr, RadialGrid::BA_PARAM_XI);
SetBounceAveragePXi(BA_xi2OverB_f1, FLUXGRIDTYPE_P1,
                   RadialGrid::BA_FUNC_XI_SQUARED_OVER_B, nullptr, 
                   RadialGrid::BA_PARAM_XI_SQUARED_OVER_B);
SetBounceAveragePXi(BA_B3_f1, FLUXGRIDTYPE_P1,
                   RadialGrid::BA_FUNC_B_CUBED, nullptr, 
                   RadialGrid::BA_PARAM_B_CUBED);
```

These are used directly in physics terms without re-computation.

---

## 8. Complete Integration Flow

### 8.1 Initialization Sequence

1. **Grid construction** (`Grid.cpp`):
   ```cpp
   Grid::Grid(RadialGrid *rg, MomentumGrid *mg, ...) {
       rgrid = rg;
       momentumGrids = new MomentumGrid*[rgrid->GetNr()];
       
       // Create BounceAverager
       bounceAverager = new BounceAverager(this, FSA, ntheta_interp_trapped, qm_trapped);
   }
   ```

2. **BounceAverager initialization** (`BounceAverager.cpp`):
   ```cpp
   BounceAverager::BounceAverager(Grid *g, FluxSurfaceAverager* fsa, ...) {
       // Set up quadrature rules
       InitializeQuadrature(q_method_trapped);
       
       // Create geometric quantity interpolators
       BOverBmin = new BounceSurfaceQuantity(grid, fsa->GetBOverBmin());
       ROverR0 = new BounceSurfaceQuantity(grid, fsa->GetROverR0());
       NablaR2 = new BounceSurfaceQuantity(grid, fsa->GetNablaR2());
       Metric = new BounceSurfaceMetric(grid, fsa->GetJacobian(), ...);
       
       // Allocate GSL workspaces
       gsl_fsolver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
       gsl_adaptive = gsl_integration_workspace_alloc(1000);
       qaws_table = gsl_integration_qaws_table_alloc(-0.5, -0.5, 0, 0);
   }
   ```

3. **Rebuild after magnetic field changes**:
   ```cpp
   void Grid::RebuildJacobians() {
       rgrid->RebuildJacobians();
       bounceAverager->Rebuild();
       RebuildBounceAveragedQuantities();
   }
   ```

### 8.2 Rebuild Process

The `Rebuild()` method updates all bounce-related quantities:

```cpp
void BounceAverager::Rebuild() {
    // Update grid resolution
    UpdateGridResolution();
    
    // Deallocate old data
    if (isTrapped != nullptr) {
        Metric->DeallocateData();
        BOverBmin->DeallocateData();
        ROverR0->DeallocateData();
        NablaR2->DeallocateData();
    }
    
    // Initialize trapped status and bounce points
    bool hasTrapped = InitializeBounceIntegralQuantities();
    
    // Determine storage requirements
    bool storeTrapped = hasTrapped && !integrateTrappedAdaptive;
    bool storePassing = !integratePassingAdaptive;
    
    // Allocate memory
    if (storePassing || storeTrapped)
        Metric->AllocateData();
    if (storeTrapped) {
        BOverBmin->AllocateData();
        ROverR0->AllocateData();
        NablaR2->AllocateData();
    }
    
    // Store data on quadrature grids
    if (storePassing)
        Metric->SetDataForPassing(ntheta_interp_passing, theta_passing);
    if (storeTrapped) {
        BOverBmin->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
        ROverR0->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
        NablaR2->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
        Metric->SetDataForTrapped(ntheta_interp_trapped, theta_trapped_ref);
    }
    
    // Calculate and store bounce-averaged metric Vp
    real_t **Vp, **Vp_fr, **Vp_f1, **Vp_f2, **VpOverP2AtZero;
    SetVpsPXi(Vp, Vp_f1, VpOverP2AtZero);
    SetVp(Vp_fr, FLUXGRIDTYPE_RADIAL);
    SetVp(Vp_f2, FLUXGRIDTYPE_P2);
    
    grid->SetVp(Vp, Vp_fr, Vp_f1, Vp_f2, VpOverP2AtZero);
}
```

### 8.3 Step-by-Step Bounce Average Calculation

For a given grid point (ir, i, j) and function F:

**Step 1**: Determine particle properties
```cpp
real_t xi0 = GetXi0(ir, i, j, fluxGridType);
bool isTrapped = BounceSurfaceQuantity::IsTrapped(ir, i, j, fluxGridType, grid);
real_t Bmin = fluxSurfaceAverager->GetBmin(ir, fluxGridType);
real_t Bmax = fluxSurfaceAverager->GetBmax(ir, fluxGridType);
```

**Step 2**: Find bounce points (if trapped)
```cpp
if (isTrapped) {
    FindBouncePoints(ir, Bmin, theta_Bmin, theta_Bmax, 
                    FSA, xi0, fluxGridType,
                    &theta_b1, &theta_b2, gsl_fsolver, geometryIsSymmetric);
} else {
    theta_b1 = 0;
    theta_b2 = 2*M_PI;
}
```

**Step 3**: Check for special cases
```cpp
// Negative-pitch trapped
if (isTrapped && xi0 < 100*realeps)
    return 0;

// Cell contains singularity
if (shouldCellAverageBounceIntegral(ir, xi_f1, xi_f2, fluxGridType))
    return EvaluateCellAveragedBounceIntegralOverP2(...);
```

**Step 4**: Choose integration method
```cpp
if (adaptive_integration) {
    // Use GSL adaptive quadrature
    if (isTrapped && integrand_is_singular)
        gsl_integration_qaws(...);
    else
        gsl_integration_qag(...);
} else {
    // Use fixed quadrature with pre-computed values
    for (it = 0; it < ntheta; it++) {
        BounceIntegral += weights[it] * Metric[it] * Function[it];
    }
}
```

**Step 5**: Normalize by phase-space volume
```cpp
real_t Vp = GetVp(ir, i, j, fluxGridType);
return BounceIntegral / Vp;
```

### 8.4 Usage in Physics Terms

Example: Electric field acceleration term

```cpp
// From documentation, Eq. (12-13)
{A_E^p} = -e {E_∥ ξ}

// Implementation would call:
real_t E_parallel_xi_avg = grid->CalculateBounceAverage(
    ir, i, j, FLUXGRIDTYPE_DISTRIBUTION,
    [](real_t xiOverXi0, real_t BOverBmin, real_t ROverR0, 
       real_t NablaR2, void *par) {
        // E_∥ is independent of position for uniform field
        // ξ = xiOverXi0 * xi0
        return E_parallel * xiOverXi0 * xi0;
    },
    &parameters, Flist
);
```

---

## 9. Key Equations

### 8.1 Fundamental Definitions

**Phase-space coordinates**:
- r: flux-surface label
- p: normalized momentum (p = γv/c)
- ξ₀: pitch at B_min: ξ₀ = (B·p)/(Bp)|_{B=B_min}

**Local pitch**:
$$\xi = \text{sgn}(\xi_0)\sqrt{1 - (1-\xi_0^2)\frac{B(\theta)}{B_{\min}}}$$

**Trapped-passing boundary**:
$$\xi_T = \sqrt{1 - \frac{B_{\min}}{B_{\max}}}$$

### 8.2 Metrics and Jacobians

**Spatial Jacobian**:
$$\mathcal{J} = \frac{1}{|\nabla\phi \cdot (\nabla\theta \times \nabla r)|}$$

**Phase-space metric**:
$$\sqrt{g} = p^2 \frac{B}{B_{\min}} \frac{\xi_0}{\xi} \mathcal{J}$$

**Phase-space volume element**:
$$\mathcal{V}' = \int_0^{2\pi} d\zeta \int_0^{2\pi} d\phi \oint d\theta \sqrt{g}$$

**Spatial volume element**:
$$V' = \int_0^{2\pi} d\phi \int_{-\pi}^{\pi} d\theta \mathcal{J}$$

### 8.3 Bounce Average Formulas

**General bounce average**:
$$\{X\} = \frac{1}{\mathcal{V}'} \int_0^{2\pi} d\zeta \int_0^{2\pi} d\phi \oint d\theta \sqrt{g} X$$

**Flux surface average**:
$$\langle Y \rangle = \frac{1}{V'} \int_0^{2\pi} d\phi \int_{-\pi}^{\pi} d\theta \mathcal{J} Y$$

**Relation between bounce and flux averages** (for passing particles):
$$\{X\} = \frac{V'}{\mathcal{V}'} \langle X \rangle$$

### 8.4 Common Bounce-Averaged Quantities

**Electric field term**:
$$\{E_\parallel \xi\} = 2\pi p^2 \xi_0 \frac{V'}{\mathcal{V}'} \frac{\langle \boldsymbol{E}\cdot\boldsymbol{B} \rangle}{B_{\min}} \times \begin{cases} 0 & \text{trapped} \\ 1 & \text{passing} \end{cases}$$

**Collision operator diffusion**:
$$\{D_C^{\xi_0\xi_0}\} = (1-\xi_0^2) \frac{\nu_D}{2} \left\{\frac{B_{\min}}{B} \frac{\xi^2}{\xi_0^2}\right\}$$

**Synchrotron radiation**:
$$\{A_S^p\} = \frac{1}{\tau_{S,\min}} p \gamma (1-\xi_0^2) \left\{\frac{B^3}{B_{\min}^3}\right\}$$

$$\{A_S^{\xi_0}\} = \frac{1}{\tau_{S,\min}} \frac{1}{\gamma} \xi_0 (1-\xi_0^2) \left\{\frac{B^2}{B_{\min}^2} \frac{\xi^2}{\xi_0^2}\right\}$$

### 8.5 Avalanche Source Term

**Rosenbluth-Putvinski source**:
$$\{C_{\text{ava}}\} = \langle n_{\text{re}} \rangle n_{\text{tot}} c r_0^2 \frac{\{B \delta(\xi-\xi^\star)\}}{\langle B \rangle} \frac{1}{p^2} \frac{\partial}{\partial p} \frac{1}{1-\gamma}$$

where:
$$\xi^\star = \sqrt{\frac{\gamma-1}{\gamma+1}}$$

The delta function bounce average requires special treatment:
$$\{B \delta(\xi-\xi^\star)\} = \text{evaluated using avalancheDeltaHat pre-computed tables}$$

---

## 10. Practical Considerations

### 10.1 Choosing Quadrature Parameters

**Number of quadrature points**:
- Default: `ntheta_interp_trapped = 30` (set in `RadialGrid.cpp`)
- Increase for higher accuracy or complex magnetic geometries
- Decrease for faster simulations with acceptable accuracy loss

**Adaptive vs Fixed**:
- **Adaptive**: Higher accuracy, slower, essential for complex integrands
- **Fixed Chebyshev**: Good balance for most cases, handles singularities well
- **Fixed Legendre**: Fastest, suitable for smooth passing particle integrals

### 10.2 Performance Tips

1. **Use PXi grids when possible**: Exploits p² scaling for 10-100× speedup
2. **Pre-compute common bounce averages**: Avoid redundant calculations
3. **Choose appropriate quadrature**: Chebyshev for trapped, Legendre for passing
4. **Cache geometric quantities**: Interpolate once, use many times
5. **Exploit symmetry**: Up-down symmetric fields reduce computation by 2×

### 10.3 Debugging Common Issues

**Problem**: Bounce points not found
- **Cause**: Magnetic field has multiple minima/maxima
- **Solution**: Check `HasMagneticFieldMultipleOptima()` in radial grid

**Problem**: NaN in bounce averages
- **Cause**: Division by zero at ξ₀ = 0 or p = 0
- **Solution**: Verify special case handling is active

**Problem**: Poor convergence near trapped-passing boundary
- **Cause**: Insufficient resolution or wrong quadrature
- **Solution**: Increase `ntheta_interp_trapped` or use adaptive quadrature

**Problem**: Asymmetric results in symmetric geometry
- **Cause**: Incorrect bounce point finding
- **Solution**: Check `geometryIsSymmetric` flag and mirror symmetry logic

### 10.4 Validation Checks

The code includes several validation checks:

```cpp
// Ensure magnetic field has single minimum and maximum
if (g->GetMomentumGrid(0)->GetNCells() > 1) {
    if (g->GetRadialGrid()->HasMagneticFieldMultipleOptima())
        FVMException("Magnetic field has multiple minima or maxima. "
                    "Incompatible with kinetic simulations.");
}

// Verify bounce points are found
if (status != GSL_SUCCESS)
    throw FVMException("Unable to find valid theta root.");

// Check for zero phase-space volume
if (!Vp)
    return 0;  // Avoid division by zero
```

---

## 11. Summary

The DREAM bounce average implementation combines:

1. **Rigorous theoretical foundation** based on drift-kinetic theory
2. **Robust numerical methods** using GSL integration and root-finding
3. **Efficient optimizations** exploiting problem structure (p² scaling, pre-computation)
4. **Careful singularity handling** for trapped-passing boundary and special points
5. **Flexible architecture** supporting multiple quadrature methods and grid types

Key strengths:
- Handles both trapped and passing particles correctly
- Robust treatment of integrable singularities
- Efficient for large-scale simulations through optimization
- Validated against analytical limits and other codes

This implementation enables DREAM to efficiently simulate runaway electron dynamics in realistic tokamak geometries while maintaining physical accuracy.

---

## References

1. DREAM paper: "DREAM: a fluid-kinetic framework for tokamak disruption runaway electron simulations"
2. GSL Library Documentation: https://www.gnu.org/software/gsl/doc/html/
3. Hirvijoki et al., Computer Physics Communications (2018)
4. Theory notes: `/data/zhzhou/DREAM/doc/notes/theory.tex`

---

## File Locations

- **Main implementation**: `/data/zhzhou/DREAM/fvm/Grid/Grid.cpp`
- **Bounce averager**: `/data/zhzhou/DREAM/fvm/Grid/BounceAverager.cpp`
- **Flux surface averager**: `/data/zhzhou/DREAM/fvm/Grid/FluxSurfaceAverager.cpp`
- **Headers**: `/data/zhzhou/DREAM/include/FVM/Grid/`
- **Theory**: `/data/zhzhou/DREAM/doc/notes/theory.tex`
- **Paper**: `/data/zhzhou/DREAM/doc/MinerU_markdown_202604281659693.md`

---

*Document generated: Based on DREAM codebase analysis and theoretical paper*
*Last updated: 2026-04-29*
