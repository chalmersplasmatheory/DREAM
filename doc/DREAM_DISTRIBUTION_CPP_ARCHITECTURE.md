# DREAM分布函数方程构建器 - distribution.cpp 架构详解

## 📋 文件概览

**文件路径**: `DREAM/src/Settings/Equations/distribution.cpp`

**核心职责**: 构建分布函数（f_hot, f_re）的动力学方程，将各种物理效应模块化地组装成FVM算子。

**设计理念**: 
- **模块化**: 每个物理效应都是独立的 `EquationTerm` 对象
- **可配置**: 通过Settings系统动态启用/禁用各项
- **可扩展**: 新增物理项只需添加新的 `AddTerm()` 调用
- **统一接口**: 所有项都遵循 FVM 框架的标准接口

---

## 🏗️ 程序结构总览

```
distribution.cpp
│
├── 1. DefineOptions_f_general()        ← 定义配置选项
│   ├── 边界条件设置
│   ├── 插值方法设置
│   ├── 物理效应开关
│   │   ├── ripple mode
│   │   ├── synchrotron mode
│   │   ├── time-varying B mode
│   │   └── quasilinear diffusion mode ⭐
│   ├── 初始分布设置
│   └── 输运模型设置
│
├── 2. ConstructEquation_f_general()    ← 构建方程（核心函数）⭐⭐⭐
│   ├── Step 1: Transient Term          ← ∂f/∂t
│   ├── Step 2: Electric Field Term     ← E·∂f/∂p
│   ├── Step 3: Pitch Scattering        ← 碰撞pitch散射
│   ├── Step 4: Ripple Scattering       ← 磁ripple（可选）
│   ├── Step 5: Time-varying B          ← dB/dt效应（可选）
│   ├── Step 6: Trapping BC             ← 捕获-通行边界
│   ├── Step 7: Slowing Down            ← 碰撞能量损失（必有）
│   ├── Step 8: Energy Diffusion        ← 碰撞能量扩散（必有）
│   ├── Step 9: Synchrotron             ← 同步辐射（可选）
│   ├── Step 10: Quasi-linear Diffusion ← 准线性扩散（可选）⭐
│   │   └── ConstructEquation_f_quasilinear() ← 封装的辅助函数
│   ├── Step 11: Radial Transport       ← 径向输运（可选）
│   ├── Step 12: External BC            ← p_max边界条件
│   ├── Step 13: Interpolation Setup    ← 通量限制器
│   ├── Step 14: Register Operator      ← 注册到方程系统
│   └── Step 15: Set Initial Value      ← 设置初始分布
│       ├── Numerical init (from file)
│       └── Maxwellian init (n0, T0)
│
├── 3. ConstructEquation_f_ripple()     ← 辅助：构建ripple项
├── 4. ConstructEquation_f_timevaryingb() ← 辅助：构建时变B项
├── 5. ConstructEquation_f_prescribed() ← 辅助：预设分布（不演化）
└── 6. ConstructEquation_f_maxwellian() ← 辅助：构造Maxwellian初值

**Note**: Quasilinear diffusion construction is now in a separate file:
- `src/Equations/Kinetic/QuasilinearDiffusionBuilder.cpp`
- Header: `include/DREAM/Equations/Kinetic/QuasilinearDiffusionBuilder.hpp`
- Function: `ConstructQuasilinearDiffusionTerm()`
```

---

## 🔧 详细函数解析

### 1️⃣ `DefineOptions_f_general()` 

**签名**:
```cpp
void SimulationGenerator::DefineOptions_f_general(
    Settings *s,           // Settings对象
    const string& mod      // 模块名称（如 "eqsys/f_hot"）
)
```

**功能**: 在Settings系统中注册所有可配置选项

**注册的选项分类**:

#### A. 数值方法配置 (Lines 48-57)

| 选项名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `boundarycondition` | int | `BC_PHI_CONST` | p_max处的边界条件类型 |
| `adv_interp/r` | int | `AD_INTERP_CENTRED` | 径向通量插值方法 |
| `adv_interp/p1` | int | `AD_INTERP_CENTRED` | 动量方向通量插值 |
| `adv_interp/p2` | int | `AD_INTERP_CENTRED` | pitch方向通量插值 |
| `adv_interp/fluxlimiterdamping` | real | 1.0 | 通量限制器阻尼系数 |

**作用**: 控制FVM离散化的数值稳定性

#### B. 物理效应开关 (Lines 59-75)

```cpp
// 磁Ripple散射
s->DefineSetting(mod + "/ripplemode", "...", 
    (int_t)OptionConstants::EQTERM_RIPPLE_MODE_NEGLECT);

// 同步辐射
s->DefineSetting(mod + "/synchrotronmode", "...", 
    (int_t)OptionConstants::EQTERM_SYNCHROTRON_MODE_NEGLECT);

// 时变磁场
s->DefineSetting(mod + "/timevaryingbmode", "...", 
    (int_t)OptionConstants::EQTERM_TIMEVARYINGB_MODE_NEGLECT);

// ⭐ 准线性扩散（我们的重点）
s->DefineSetting(mod + "/quasilinearmode", "...", 
    (int_t)OptionConstants::QL_DIFFUSION_MODE_NEGLECT);
```

**准线性扩散的详细配置**:

```cpp
// 是否使用预计算矩阵
s->DefineSetting(mod + "/quasilinear/use_precomputed_matrix", "...", (int_t)0);

// HDF5文件路径
s->DefineSetting(mod + "/quasilinear/precomputed_file", "...", std::string(""));

// 波振幅（归一化到 n_e m_e c^2）
s->DefineSetting(mod + "/quasilinear/amplitude", "...", (real_t)1e-10);

// 波谱类型
s->DefineSetting(mod + "/quasilinear/spectrum_type", "...", 
    (int_t)OptionConstants::WAVE_SPECTRUM_UNIFORM);

// k网格分辨率
s->DefineSetting(mod + "/quasilinear/num_k", "...", (int_t)100);
s->DefineSetting(mod + "/quasilinear/num_ktheta", "...", (int_t)20);

// k范围
s->DefineSetting(mod + "/quasilinear/k_min", "...", (real_t)35.0);
s->DefineSetting(mod + "/quasilinear/k_max", "...", (real_t)45.0);
s->DefineSetting(mod + "/quasilinear/ktheta_min", "...", (real_t)0.1);
s->DefineSetting(mod + "/quasilinear/ktheta_max", "...", (real_t)0.3);

// 谐波模式
s->DefineSetting(mod + "/quasilinear/harmonic_mode", "...", 
    (int_t)OptionConstants::QL_HARMONIC_N_MINUS_1);
```

#### C. 初始分布配置 (Lines 79-82)

```cpp
DefineDataR(mod, s, "n0");      // 密度剖面 n(r)
DefineDataR(mod, s, "T0");      // 温度剖面 T(r)
DefineDataR2P(mod, s, "init");  // 数值指定的 f(r,p,ξ)
```

#### D. 其他配置 (Lines 84-91)

```cpp
DefineOptions_Transport(mod, s, true);  // 径向输运模型
DefineDataTR2P(mod, s, "f_prescribed"); // 预设分布（时间依赖）
s->DefineSetting(mod + "/fullIonJacobian", "...", (bool)true); // Jacobian精度
```

---

### 2️⃣ `ConstructEquation_f_general()` ⭐⭐⭐ 核心函数

**签名**:
```cpp
FVM::Operator *SimulationGenerator::ConstructEquation_f_general(
    Settings *s,                    // Settings对象
    const string& mod,              // 模块名称
    EquationSystem *eqsys,          // 方程系统
    len_t id_f,                     // 分布函数的ID
    FVM::Grid *grid,                // 动量网格
    enum OptionConstants::momentumgrid_type gridtype, // 网格类型
    CollisionQuantityHandler *cqty, // 碰撞量处理器
    bool addExternalBC,             // 是否添加外部边界条件
    bool addInternalBC,             // 是否添加内部边界条件
    FVM::Operator **transport,      // 输运算子指针
    TransportAdvectiveBC **advective_bc,   // 对流边界条件
    TransportDiffusiveBC **diffusive_bc,   // 扩散边界条件
    RipplePitchScattering **ripple_Dxx,    // Ripple项指针
    SynchrotronTerm **synchrotron,         // 同步辐射项指针
    TimeVaryingBTerm **timevaryingb,       // 时变B项指针
    bool rescaleMaxwellian          // 是否重缩放Maxwellian
)
```

**返回值**: `FVM::Operator*` - 构建好的动力学算子

**完整流程**:

#### Step 1: 创建算子对象 (Line 121)
```cpp
FVM::Operator *eqn = new FVM::Operator(grid);
```
- 创建一个空的FVM算子
- 关联到指定的动量网格

#### Step 2: 添加瞬态项 (Line 124)
```cpp
eqn->AddTerm(new FVM::TransientTerm(grid, id_f));
```
- 对应 ∂f/∂t
- **所有方程必须包含此项**

#### Step 3: 电场加速项 (Lines 126-150)

根据是否是简化方程选择不同的实现：

```cpp
bool isReducedEquation = 
    (gridtype == OptionConstants::MOMENTUMGRID_TYPE_PXI &&
     grid->GetMomentumGrid(0)->GetNp2() == 1);

if (isReducedEquation) {
    // 简化方程：用扩散形式表示电场
    eqn->AddTerm(new ElectricFieldDiffusionTerm(
        grid, cqty, eqsys->GetUnknownHandler(), withFullIonJacobian
    ));
} else {
    // 完整3D方程：用平流形式表示电场
    eqn->AddTerm(new ElectricFieldTerm(
        grid, eqsys->GetUnknownHandler(), gridtype
    ));
}
```

**物理意义**: E·∂f/∂p（电场加速电子）

#### Step 4: Pitch散射项 (Lines 153-157)
```cpp
eqn->AddTerm(new PitchScatterTerm(
    grid, cqty, gridtype,
    eqsys->GetUnknownHandler(),
    withFullIonJacobian
));
```

**物理意义**: ν_D · ∂/∂ξ[(1-ξ²)∂f/∂ξ]（碰撞导致的pitch角度散射）

#### Step 5: 磁Ripple散射 (Lines 160-161)
```cpp
if ((*ripple_Dxx = ConstructEquation_f_ripple(s, mod, grid, gridtype)) != nullptr)
    eqn->AddTerm(*ripple_Dxx);
```

**特点**:
- 可选，由 `ripplemode` 控制
- 托卡马克波纹场导致的额外pitch散射
- 通过辅助函数 `ConstructEquation_f_ripple()` 构建

#### Step 6: 时变磁场效应 (Lines 164-165)
```cpp
if ((*timevaryingb = ConstructEquation_f_timevaryingb(s, mod, grid)) != nullptr)
    eqn->AddTerm(*timevaryingb);
```

**特点**:
- 可选，由 `timevaryingbmode` 控制
- dB/dt ≠ 0 时的绝热压缩效应
- 通过辅助函数 `ConstructEquation_f_timevaryingb()` 构建

#### Step 7: 捕获-通行边界条件 (Lines 169-170)
```cpp
if(grid->HasTrapped())
    eqn->AddBoundaryCondition(new FVM::BC::PXiInternalTrapping(grid, eqn));
```

**物理意义**: 处理 trapped/passing 粒子在非均匀磁场中的转换

#### Step 8: Slowing Down项 (Lines 175-179) ⚡ **始终存在**
```cpp
eqn->AddTerm(new SlowingDownTerm(
    grid, cqty, gridtype, 
    eqsys->GetUnknownHandler(),
    withFullIonJacobian
));
```

**物理意义**: -∂/∂p[ν_s · p · f]（碰撞导致的能量损失）

**重要性**: 这是Fokker-Planck方程的核心项之一，描述高能电子通过与背景等离子体碰撞而减速

#### Step 9: 能量扩散项 (Lines 182-186) ⚡ **始终存在**
```cpp
eqn->AddTerm(new EnergyDiffusionTerm(
    grid, cqty, gridtype,
    eqsys->GetUnknownHandler(),
    withFullIonJacobian
));
```

**物理意义**: ∂/∂p[D_pp · ∂f/∂p]（碰撞导致的能量扩散）

**重要性**: 与Slowing Down项一起构成完整的碰撞算子

#### Step 10: Synchrotron (Lines 189-194)
```cpp
enum OptionConstants::eqterm_synchrotron_mode synchmode =
    (enum OptionConstants::eqterm_synchrotron_mode)s->GetInteger(mod + "/synchrotronmode");

if (synchmode == OptionConstants::EQTERM_SYNCHROTRON_MODE_INCLUDE) {
    *synchrotron = new SynchrotronTerm(grid, gridtype);
    eqn->AddTerm(*synchrotron);
}
```

**特点**:
- 可选，由 `synchrotronmode` 控制
- 同步辐射导致的能量损失和pitch变化
- 对高能电子特别重要

#### Step 11: 准线性扩散 ⭐⭐⭐ (Lines 196-201)

这是**我们实现的核心部分**，现在已封装为独立的Builder函数，并使用与其他项一致的指针参数模式！

```cpp
// Quasilinear diffusion from external waves
if (quasilinear != nullptr) {
    *quasilinear = ConstructQuasilinearDiffusionTerm(s, mod, grid);
    if (*quasilinear != nullptr)
        eqn->AddTerm(*quasilinear);
}
```

**关键改进**:
- ✅ **独立文件**: 所有准线性扩散构建逻辑在 `src/Settings/Equations/Kinetic/QuasilinearDiffusionBuilder.cpp`
- ✅ **使用 new 模式**: 与 ripple/synchrotron/timevaryingb 保持一致，通过指针参数返回
- ✅ **外部可访问**: 指针存储在 `OtherQuantityHandler::eqn_terms` 中，可用于诊断和输出
- ✅ **两种模式支持**: 预计算矩阵模式和实时计算模式都在Builder函数中处理
- ✅ **易于维护**: 与distribution.cpp解耦，便于单独测试和修改
- ✅ **清晰的架构**: Builder模式 + 指针参数模式，职责分离

**Builder函数位置**: 
- Header: [QuasilinearDiffusionBuilder.hpp](file:///data/zhzhou/DREAM/include/DREAM/Equations/Kinetic/QuasilinearDiffusionBuilder.hpp)
- Implementation: [QuasilinearDiffusionBuilder.cpp](file:///data/zhzhou/DREAM/src/Equations/Kinetic/QuasilinearDiffusionBuilder.cpp)

**外部访问示例** (在 OtherQuantityHandler 中):
```cpp
// 可以访问准线性扩散项来计算诊断量
if (tracked_terms->f_hot_quasilinear != nullptr) {
    // 计算准线性加热功率、扩散系数等
    auto D_pp = tracked_terms->f_hot_quasilinear->getD_pp(ir, i, j);
    // ...
}
```

**详细实现见下方第5节**

#### Step 12: 径向输运 (Lines 312-327)
```cpp
bool hasTransport = ConstructTransportTerm(
    eqn, mod, grid, gridtype, eqsys, s, true, false, 
    advective_bc, diffusive_bc
);

if (hasTransport && transport != nullptr) {
    *transport = new FVM::Operator(grid);
    ConstructTransportTerm(
        *transport, mod, grid, gridtype, eqsys, s, true, false
    );
}
```

**特点**:
- Rechester-Rosenbluth输运等
- 描述湍流导致的径向扩散
- 可以单独输出为一个算子

#### Step 13: 外部边界条件 (Lines 331-339)
```cpp
if (addExternalBC) {
    enum FVM::BC::PXiExternalLoss::bc_type bc =
        (enum FVM::BC::PXiExternalLoss::bc_type)s->GetInteger(mod + "/boundarycondition");
    
    eqn->AddBoundaryCondition(new FVM::BC::PXiExternalLoss(
        grid, eqn, id_f, nullptr,
        FVM::BC::PXiExternalLoss::BOUNDARY_KINETIC, bc
    ));
}
```

**物理意义**: p = p_max 处的粒子逃逸边界条件

#### Step 14: 设置插值方法 (Lines 342-366)
```cpp
// 读取插值配置
enum FVM::AdvectionInterpolationCoefficient::adv_interpolation adv_interp_r =
    (enum ...)s->GetInteger(mod + "/adv_interp/r");
// ... 同样读取 p1, p2 方向

// 设置通量限制器
eqn->SetAdvectionInterpolationMethod(
    adv_interp_r, adv_jac_mode_r, FVM::FLUXGRIDTYPE_RADIAL, 
    id_f, fluxLimiterDamping
);
eqn->SetAdvectionInterpolationMethod(
    adv_interp_p1, adv_jac_mode_p1, FVM::FLUXGRIDTYPE_P1, 
    id_f, fluxLimiterDamping
);
eqn->SetAdvectionInterpolationMethod(
    adv_interp_p2, adv_jac_mode_p2, FVM::FLUXGRIDTYPE_P2, 
    id_f, fluxLimiterDamping
);
```

**作用**: 防止数值振荡，保证稳定性

#### Step 15: 内部边界条件 (Lines 370-371)
```cpp
if (addInternalBC)
    eqn->SetAdvectionBoundaryConditions(
        FVM::FLUXGRIDTYPE_P1, 
        FVM::AdvectionInterpolationCoefficient::AD_BC_MIRRORED, 
        FVM::AdvectionInterpolationCoefficient::AD_BC_DIRICHLET
    );
```

**物理意义**: p = 0 处的镜像边界条件

#### Step 16: 注册算子 (Line 373)
```cpp
eqsys->SetOperator(id_f, id_f, eqn, desc);
```

**关键步骤**:
- 将构建好的算子注册到方程系统
- 后续时间积分会使用这个算子
- `desc` 是描述字符串（"3D kinetic equation" 或 "Reduced kinetic equation"）

#### Step 17: 设置初始值 (Lines 376-411)

有两种方式设置初始分布：

**方式A: 数值指定** (Lines 380-388)
```cpp
len_t nx[3];
if (s->GetRealArray(mod + "/init/x", 3, nx, false) != nullptr) {
    // 从输入文件加载数值分布
    FVM::Interpolator3D *interp = LoadDataR2P(mod, s, "init");
    enum FVM::Interpolator3D::momentumgrid_type momtype = 
        GetInterp3DMomentumGridType(gridtype);
    const real_t *init = interp->Eval(grid, momtype);
    
    eqsys->SetInitialValue(id_f, init);
    
    delete [] init;
    delete interp;
}
```

**方式B: Maxwellian分布** (Lines 389-411)
```cpp
else {
    real_t *n0 = LoadDataR(mod, grid->GetRadialGrid(), s, "n0");
    real_t *T0 = LoadDataR(mod, grid->GetRadialGrid(), s, "T0");
    
    ConstructEquation_f_maxwellian(id_f, eqsys, grid, n0, T0, rescaleMaxwellian);
    
    delete [] T0;
    delete [] n0;
}
```

通过辅助函数 `ConstructEquation_f_maxwellian()` 构造：
```cpp
f(p) = n₀ / (4π m_e³ c³ θ³) · exp(-(γ-1)/θ)
```
其中 θ = T/(m_e c²)

**可选的重缩放** (Lines 536-540):
```cpp
if(rescaleMaxwellian){
    real_t normalizationFactor = n0[ir]/grid->IntegralMomentumAtRadius(ir,f);
    for (len_t i = 0; i < np1*np2; i++)
        f[i] *= normalizationFactor;
}
```
确保数值积分精确等于 n₀

---

### 3️⃣ 辅助函数

#### `ConstructEquation_f_ripple()` (Lines 420-457)

**功能**: 构建磁ripple散射项

**流程**:
1. 检查 `ripplemode` 是否启用
2. 读取线圈数 `ncoils` 或波纹幅度 `deltaCoils`
3. 加载傅里叶模式 (m, n) 和 dB/B 剖面
4. 创建 `RipplePitchScattering` 对象

#### `ConstructEquation_f_timevaryingb()` (Lines 462-475)

**功能**: 构建时变磁场项

**流程**:
1. 检查 `timevaryingbmode` 是否启用
2. 加载 d(ln B₀)/dt 剖面
3. 创建 `TimeVaryingBTerm` 对象

#### `ConstructQuasilinearDiffusionTerm()` ⭐ (Separate File)

**文件位置**: 
- Header: `include/DREAM/Equations/Kinetic/QuasilinearDiffusionBuilder.hpp`
- Implementation: `src/Equations/Kinetic/QuasilinearDiffusionBuilder.cpp`

**功能**: 构建准线性扩散项（独立Builder函数）

**签名**:
```cpp
QuasilinearDiffusionTerm *DREAM::ConstructQuasilinearDiffusionTerm(
    Settings *s, 
    const std::string& mod, 
    FVM::Grid *grid
)
```

**返回值**: 
- `QuasilinearDiffusionTerm*` - 如果启用了准线性扩散
- `nullptr` - 如果未启用 (`quasilinearmode == QL_DIFFUSION_MODE_NEGLECT`)

**详细流程**: See implementation in [QuasilinearDiffusionBuilder.cpp](file:///data/zhzhou/DREAM/src/Equations/Kinetic/QuasilinearDiffusionBuilder.cpp)

The function supports two modes:
1. **Pre-computed matrix mode**: Loads diffusion coefficients from HDF5 file
2. **On-the-fly computation mode**: Calculates coefficients using wave-particle interaction theory

**设计优势**:

| 特性 | 说明 |
|------|------|
| **模块化** | 独立的文件，与distribution.cpp解耦 |
| **可测试性** | 可以单独编译和测试Builder逻辑 |
| **灵活性** | 支持两种模式（预计算和实时计算），通过Settings切换 |
| **错误处理** | 在预计算模式下验证HDF5文件路径，抛出有意义的异常 |
| **日志输出** | 提供详细的初始化信息，方便调试 |
| **架构清晰** | Builder模式，职责分离 |

**调用方式**:
```cpp
// 在 ConstructEquation_f_general() 中
#include "DREAM/Equations/Kinetic/QuasilinearDiffusionBuilder.hpp"

QuasilinearDiffusionTerm *qlTerm;
if ((qlTerm = ConstructQuasilinearDiffusionTerm(s, mod, grid)) != nullptr)
    eqn->AddTerm(qlTerm);  // ← 添加到方程
```

This follows a clean separation of concerns pattern!

#### `ConstructEquation_f_prescribed()` (Lines 480-492)

**功能**: 构建预设分布方程（不演化）

**用途**: 当分布函数被外部指定且不随时间演化时使用

#### `ConstructEquation_f_maxwellian()` (Lines 502-546)

**功能**: 构造Maxwellian初始分布

**详细流程**:
```cpp
for each radial point ir:
    1. 获取动量网格点 p[i,j]
    2. 计算相对论Maxwellian:
       f[i,j] = Constants::RelativisticMaxwellian(p, n0[ir], T0[ir])
    3. 检查网格分辨率（f(p=0)不应为0）
    4. 可选：重缩放使积分精确等于n0
```

**警告机制** (Lines 526-533):
```cpp
if (f[0] == 0)
    DREAM::IO::PrintWarning(
        DREAM::IO::WARNING_P_UNDERRESOLVED,
        "The momentum grid resolution is so low that the initial "
        "Maxwellian cannot be accurately represented on the grid."
    );
```

---

## 🎯 设计模式分析

### 1. 组合模式 (Composite Pattern)

```
FVM::Operator
├── TransientTerm
├── ElectricFieldTerm
├── PitchScatterTerm
├── SlowingDownTerm
├── EnergyDiffusionTerm
├── SynchrotronTerm (optional)
├── QuasilinearDiffusionTerm (optional)
└── TransportTerm (optional)
```

**优势**: 
- 每个物理项独立开发、测试
- 易于添加新物理效应
- 运行时动态组合

### 2. 策略模式 (Strategy Pattern)

准线性扩散的两种实现：
- **预计算策略**: 从HDF5快速加载
- **实时计算策略**: 灵活但较慢

通过 `use_precomputed_matrix` 开关切换

### 3. 工厂模式 (Factory Pattern)

`ConstructEquation_f_general()` 是一个工厂函数，根据配置动态创建不同的项组合

---

## 📊 数据流图

```
Python脚本 (run_ql_precomputed_simulation.py)
    ↓ 设置参数
Settings对象 (DREAMSettings)
    ↓ 保存为HDF5
    ↓
C++读取Settings
    ↓
DefineOptions_f_general()  ← 注册选项
    ↓
ConstructEquation_f_general()  ← 构建方程
    ├─ 创建空Operator
    ├─ 逐项AddTerm()
    │   ├─ TransientTerm
    │   ├─ ElectricFieldTerm
    │   ├─ Collisions (Pitch + Slowing + Energy)
    │   ├─ SynchrotronTerm (可选)
    │   ├─ QuasilinearDiffusionTerm (可选) ⭐
    │   │   ├─ 预计算模式: 从HDF5加载
    │   │   └─ 实时模式: 创建spectrum+dispersion+resonance
    │   └─ Transport (可选)
    ├─ 设置边界条件
    ├─ 设置插值方法
    ├─ 注册到EquationSystem
    └─ 设置初始值 (Maxwellian或数值)
    ↓
返回 Operator*
    ↓
EquationSystem开始时间演化
    ↓
LinearImplicitSolver
    ├─ Rebuild(t, dt)  ← 重建所有项的系数
    ├─ Assemble Jacobian  ← 组装Jacobian矩阵
    └─ Solve  ← 求解线性系统
    ↓
输出结果 (f_re, n_re, ...)
```

---

## 🔑 关键要点总结

### ✅ 架构优势

1. **模块化**: 每个物理效应独立，易于维护和扩展
2. **可配置**: 通过Settings系统灵活控制
3. **统一接口**: 所有项都遵循FVM标准
4. **性能优化**: 预计算模式避免重复计算

### ⚠️ 注意事项

1. **添加顺序**: 某些项可能有依赖关系（如边界条件应在最后）
2. **内存管理**: 注意 `new` 的对象由谁负责 `delete`
3. **网格一致性**: 所有项必须使用相同的网格
4. **数值稳定性**: 合理设置通量限制器和时间步长

### 🚀 扩展指南

要添加新的物理项：

1. 创建新的 `EquationTerm` 子类
2. 在 `DefineOptions_f_general()` 中添加配置选项
3. 在 `ConstructEquation_f_general()` 中添加 `AddTerm()` 调用
4. （可选）创建辅助构造函数
5. 更新Python接口

**示例**: 添加新的辐射损失项
```cpp
// 1. 在DefineOptions中添加
s->DefineSetting(mod + "/newradiationmode", "...", (int_t)MODE_NEGLECT);

// 2. 在ConstructEquation中添加
if (radiation_mode == MODE_INCLUDE) {
    eqn->AddTerm(new NewRadiationTerm(grid, ...));
}
```

---

## 📚 相关文档

- [DREAM_QUASILINEAR_INTEGRATION_FLOW.md](file:///data/zhzhou/DREAM_QUASILINEAR_INTEGRATION_FLOW.md) - 准线性扩散集成流程详解
- [DREAM_QUASILINEAR_DIFFUSION_TECHNICAL_DETAILS.md](file:///data/zhzhou/DREAM_QUASILINEAR_DIFFUSION_TECHNICAL_DETAILS.md) - 技术细节
- FVM框架文档: `DREAM/fvm/doc/`

---

## 🎓 学习建议

1. **先理解FVM框架**: 阅读 `FVM/Equation/Operator.hpp` 了解算子结构
2. **跟踪一个简单项**: 例如 `SynchrotronTerm`，看它如何实现
3. **调试技巧**: 在 `AddTerm()` 后加断点，观察算子状态
4. **对比实现**: 比较预计算模式和实时模式的差异

这个文件是理解DREAM如何组织动力学方程的**核心入口**，掌握它对开发和调试至关重要！
