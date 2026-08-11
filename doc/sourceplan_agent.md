# AvalancheSourceDirect — 直接数值积分雪崩源项

## 前因后果：为什么需要直接数值积分

### 原有 RP 雪崩源项的根本问题

DREAM 原有的动力学雪崩源项基于 **Rosenbluth‑Putvinski (1997)** 模型，其核心假设包括：

- **入射电子完全平行于磁场**：即种子逃逸电子的投掷角 \(|\xi| = 1\)，垂直于磁场的速度分量为零。
- **入射电子能量无穷大**：以此得到解析的动量导数项 \(\partial/\partial p [1/(1-\gamma)]\)。

在这些假设下，次级电子的投掷角被严格冻结在由能量决定的单一值 \(\xi^* = \sqrt{(\gamma-1)/(\gamma+1)}\)，且源项的形状函数 \(S(r,p)\) 与种子分布的角度细节无关。

### RP 源项的完整公式

RP 模型首先定义临界电场（决定雪崩能否发生的最低电场）：

\[
E_c = \frac{4\pi e^3 n_e \ln\Lambda}{m_e c^2} \quad\text{(完全电离)},
\qquad
E_c = \frac{2\pi e^3 (n_e + n_b) \ln\Lambda}{m_e c^2} \quad\text{(包含束缚电子)}
\]

特征碰撞时间：

\[
\tau = \frac{m_0^2 c^3}{4\pi n_e e^4 \ln\Lambda} = \frac{m_0 c}{e E_c}
\]

**雪崩增长率**——RP 模型的核心结果（式 18，普适拟合公式）：

\[
\boxed{
\frac{1}{j_{\rm RA}}\frac{\partial j_{\rm RA}}{\partial t} =
\frac{1}{\tau \ln\Lambda}
\sqrt{\frac{\pi\gamma}{3(Z_{\rm eff}+5)}}
\left(\frac{E}{E_c} - 1\right)
\left[
1 - \frac{E_c}{E}
\frac{4\pi (Z_{\rm eff}+1)^2}{3\gamma (Z_{\rm eff}+5)
\left( E^2/E_c^2 + 4/\gamma^2 - 1 \right)}
\right]^{-1/2}
}
\]

其中 \(\gamma(\epsilon) \approx \left(1 + 1.46\sqrt{\epsilon} + 1.72\epsilon\right)^{-1}\) 为环径比相关的新经典因子，\(\epsilon = r/R\)。

**雪崩源项的动量空间微分形式**——这是直接与动理学方程相关的表达式：

\[
\boxed{
S_{\rm RP}(p,\xi) =
\frac{n_{\rm RE}}{4\pi\tau_c\ln\Lambda}\,
\delta(\xi - \xi_0)\,
\frac{m_e^3 c^3}{p^2}\,
\frac{\partial}{\partial p}\left(\frac{1}{1-\gamma}\right)
}
\]

其中

\[
\xi_0 = \sqrt{\frac{\gamma-1}{\gamma+1}}, \qquad
\gamma = \sqrt{1+p^2}
\]

该源项的物理含义：所有次级电子都以**固定投掷角** \(\xi_0\) 产生，角度无展宽。在 DREAM 等流体代码中，该源项通常以通量面平均形式实现：

\[
\{C_{\rm ava}\} =
\langle n_{\rm re}\rangle\,
n_{\rm tot}\,
c r_0^2\,
\frac{\{B\delta(\xi - \xi^\star)\}}{\langle B\rangle}\,
\frac{1}{p^2}\,\frac{\partial}{\partial p}\frac{1}{1-\gamma},
\qquad
\xi^\star = \sqrt{\frac{\gamma-1}{\gamma+1}}
\]

RP 模型的其他极限形式：

- **大电场极限**（忽略投掷角散射）：

\[
\frac{\partial n_r}{\partial t} \approx \frac{n_r}{2\tau \ln\Lambda} \left(\frac{E}{E_c} - 1\right)
\]

- **中等电场、强投掷角散射**（以 \(\beta\) 表示）：

\[
\frac{1}{n_r}\frac{\partial n_r}{\partial t} = \frac{1}{\tau\ln\Lambda} \left(\frac{\pi\beta}{2}\right)^{1/2}
\left(1 + \frac{8\pi}{9\beta}\right)^{1/2},
\quad
\beta = \frac{2}{3(Z_{\rm eff}+1)}\left(\frac{E}{E_c}\right)^2
\]

- **近阈值行为**（\(E\approx E_c\)）：

\[
\frac{1}{n_r}\frac{\partial n_r}{\partial t} = \frac{1}{\tau(Z_{\rm eff}+1)\ln\Lambda}
\frac{E}{E_c}\left(\frac{E}{E_c} - 1\right)
\]

**RP 模型的关键物理假定**：

1. **种子投掷角**：假设所有已存在的逃逸电子严格平行于磁场（\(|\xi|=1\)）。
2. **能量转移**：仅考虑高能种子（\(\gamma_1 \gg 1\)），忽略有限能量修正。
3. **角度分布**：次级电子的投掷角被严格固定在 \(\xi_0\)，无展宽，不考虑俘获效应。
4. **碰撞过程**：使用 Møller 截面的小角度+高能极限。

这些假设正是后续修正（如 Nilsson 模型、RE_PRL 的修正理论以及本工作开发直接数值积分方案）试图放松的。RP 模型在种子投掷角分布展宽时系统性高估雪崩产生率，因此大 pitch 角散射场景下需要更精确的源项。

### 这一假设在真实等离子体中的缺陷

在破裂过程中，逃逸电子会通过 **准线性扩散**（如哨声波、CAE 等不稳定性散射）或 **大角度碰撞** 获得非零投掷角。当种子电子具有有限投掷角时，通过 Møller 散射产生的次级电子 **也会继承展宽的投掷角分布**，其中一部分会进入 **俘获区**（\(\xi < \xi_T\)）。被磁镜捕获的电子无法持续从平行电场中获得能量，对雪崩没有贡献。

RP 源项完全忽略这种展宽，总是假设所有次级电子都产生在 \(\xi \approx 1\) 的方向上，从而 **系统性地高估了雪崩增长率**，特别是当种子分布因散射而展宽时。

### 直接数值积分方案的设计思路

`AvalancheSourceDirect` 方案正是为了解决上述问题而设计的。核心思路是：

- **放弃解析近似，直接对碰撞积分进行数值计算**。使用完整的 Møller 截面以及解析的角分布函数 \(\Pi\)（含 \(1/\sqrt{\cdots}\) 奇点），对种子分布 \(f_{\rm re}(p_1,\xi_1)\) 做二维数值积分，得到目标点 \((p,\xi)\) 的雪崩源项。
- **解析处理角度奇点**：角分布 \(\Pi\) 中的平方根奇点通过 \(\arcsin\) 变换解析积分，避免了数值发散，同时完整保留了有限投掷角导致的角展宽。
- **自然包含俘获边界**：积分区域显式限制在通行粒子区域（\(\xi > \xi_{\rm trapped}\)），俘获区内的电子不参与雪崩产生，这一筛选直接嵌入了源项计算，而不是在后期做对称平均等后处理修正。
- **无任何基函数截断误差**：与勒让德谱展开不同，直接数值积分不对角度做全局多项式展开，因此不会因俘获边界处的导数不连续或奇点附近的 Gibbs 振荡而产生较大误差，在大投掷角区域仍能保持精度。

### 两种源项对比

| 特性 | RP 源项 (原) | 直接数值积分源项 (新) |
|------|-------------|---------------------|
| 种子投掷角假设 | \(\xi = 1\) | 任意分布 |
| 次级电子角度分布 | 冻结为单一 \(\xi^*\) | 含展宽的完整分布 |
| 俘获电子处理 | 无（后处理对称平均） | 积分区域直接排除 |
| 大投掷角准确性 | 系统性高估 | 正确反映抑制 |

**一句话**：由于 RP 源项假设所有种子电子都平行于磁场，导致在大角度散射存在时严重高估雪崩电子产生率；直接数值积分法通过保留完整的碰撞运动学和角度展宽，能够自洽地捕捉种子投掷角展宽带来的雪崩抑制作用，为 DREAM 提供了更精确的雪崩源项选项。

---

## 物理模型

### 雪崩源项的数学形式

逃逸电子雪崩过程的动理学方程为：

\[
\left(\frac{\partial f_{\rm re}}{\partial t}\right)_{\rm ava}(p,\xi)
= n_{\rm tot} \iint f_{\rm re}(p_1,\xi_1) \,
K(p_1,\xi_1 \rightarrow p,\xi) \; p_1^2 dp_1 \, d\xi_1
\]

其中核函数 \(K\) 分解为能量部分和角度部分：

\[
K(p_1,\xi_1 \rightarrow p,\xi) = A(p_1, p) \cdot \Pi(\xi; p_1, \xi_1, p)
\]

### Møller 截面项 \(A(p_1, p)\)

描述初级电子 \((p_1, \xi_1)\) 通过 Møller 散射产生次级电子 \((p, \xi)\) 的能量分布。

定义能量分数：

\[
\nu = \frac{\gamma - 1}{\gamma_1 - 1}, \quad
\gamma = \sqrt{1+p^2}, \quad
\gamma_1 = \sqrt{1+p_1^2}
\]

运动学约束：\(0 < \nu < 1\)，即 \(\gamma < 2\gamma_1 - 1\)。

Møller 微分截面（相对论形式）：

\[
\frac{d\sigma}{d\nu} = \frac{\gamma_1^2}{(\gamma_1-1)^2(\gamma_1+1)}
\left[ \frac{1}{\nu^2(1-\nu)^2} - \frac{3}{\nu(1-\nu)} + \frac{(\gamma_1-1)^2}{\gamma_1^2}\left(1 + \frac{1}{\nu(1-\nu)}\right) \right]
\]

变换到 \(\frac{d\sigma}{dp}\)：

\[
\frac{d\sigma}{dp} = \frac{d\sigma}{d\nu} \cdot \frac{d\nu}{dp}, \quad
\frac{d\nu}{dp} = \frac{p}{\gamma(\gamma_1-1)}
\]

最终：

\[
A(p_1, p) = n_{\rm tot} \cdot v_e(p_1) \cdot \frac{d\sigma}{dp} \cdot \frac{2\pi}{c} \cdot r_0^2 c
\]

其中 \(r_0 = e^2/(4\pi\epsilon_0 m_e c^2)\) 为经典电子半径，\(r_0^2 c\) 为整体前置因子。

### 角度分布项 \(\Pi(\xi; p_1, \xi_1, p)\)

采用 Boozer (2015) 补充材料中的公式。

定义相对论因子：

\[
\gamma_{\rm rel} = \frac{(\gamma_1+1)(\gamma-1)}{(\gamma_1-1)(\gamma+1)},
\quad \xi_* = \sqrt{\gamma_{\rm rel}}
\]

角度分布是以 \(\xi_{1c}\) 为中心、宽度为 \(\xi_2\) 的反余弦型分布：

\[
\xi_{1c} = \xi_1 \cdot \xi_*, \quad
\xi_2 = \sqrt{1-\xi_1^2} \cdot \sqrt{1-\xi_*^2}
\]

分布函数的解析形式：

\[
\Pi(\xi; p_1, \xi_1, p) = \frac{1}{\pi \sqrt{\xi_2^2 - (\xi - \xi_{1c})^2}},
\quad |\xi - \xi_{1c}| < \xi_2
\]

在目标 ξ 单元 \([\xi_{\rm low}, \xi_{\rm high}]\) 上的解析积分为：

\[
\int_{\xi_{\rm low}}^{\xi_{\rm high}} \Pi(\xi) \, d\xi
= \frac{1}{\pi} \left[ \arcsin\!\left(\frac{\xi_{\rm high} - \xi_{1c}}{\xi_2}\right) - \arcsin\!\left(\frac{\xi_{\rm low} - \xi_{1c}}{\xi_2}\right) \right]
\]

### 相空间积分区域的限制

为减少计算量并符合物理，只对以下区域做积分：

- **动量的下限**：源电子 \(p_1 > p_{\rm cut}\)（由用户通过 `pCutAvalanche` 指定）
- **投掷角的限制**：源和目标电子均只考虑**通行粒子**（passing particles），即 \(\xi > \xi_{\rm trapped}\)，其中捕获角余弦值由磁场位形决定：

\[
\xi_{\rm trapped}(r) = \sqrt{1 - \frac{B_{\rm min}(r)}{B_{\rm max}(r)}}
\]

捕获粒子被困在磁镜中，不参与雪崩过程。

A. H. Boozer, Phys. Plasmas 22, 032504 (2015)

---

## 数值实现

### 总体架构

雪崩源项 `AvalancheSourceDirect` 继承 `FVM::EquationTerm`，是一个 **f_re 方程上的对角算子（f_re→f_re）**。

物理上，这是一个动量空间内的卷积积分。通过**预计算稀疏核矩阵**，将卷积转化为运行时的一次稀疏矩阵-向量乘法，大幅提高效率。

### 预计算流程

核矩阵 \(K_{p,\xi \leftarrow p_1,\xi_1}\) 在网格初始化时预计算，每个径向位置对应一个 PETSc 稀疏矩阵。

```
GridRebuilt()
  ├─ AllocateMatrices()           ── 获取网格尺寸 Np1, Np2, NCells
  └─ BuildKernelMatrices()
       └─ BuildKernelMatrixForRadialIndex(ir)
            ├─ Phase 1: 非零模式计数
            └─ Phase 2: 填充矩阵元素
```

#### Phase 1 — 非零模式计数

遍历所有可能的 (源 → 目标) 单元组合，筛选出有效的耦合：

1. **能量阈值**：\(\gamma < 2\gamma_1 - 1\)（Møller 运动学上界）
2. **源动量下限**：\(p_1 > p_{\rm cut}\)
3. **目标投掷角上限**：目标 ξ 单元上界 \(> \xi_{\rm trapped}\)
4. **源投掷角上限**：源 ξ 单元中心 \(> \xi_{\rm trapped}\)
5. **角度窗口重叠**：\([\xi_{1c}-\xi_2, \xi_{1c}+\xi_2]\) 与目标 ξ 单元相交
6. **相空间体积非零**：源和目标单元的 Vp > 0

每找到一个有效耦合，对应行的非零计数加一。

#### Phase 2 — 填充矩阵元素

对每个有效耦合计算矩阵元数值：

1. **Møller 截面**：`ComputeA(p1, p, 1.0)`，ntot=1 运行时缩放
2. **角度窗口参数**：`ComputeAngularSupport()` 计算 \(\xi_{1c}, \xi_2\)
3. **角度解析积分**：\([a,b] = [\xi_{\rm low}, \xi_{\rm high}] \cap [\xi_{1c}-\xi_2, \xi_{1c}+\xi_2]\)，然后用 asin 公式计算
4. **体积归一化**：考虑 Vp 权重和动量空间体积元

\[
K_{ij} = A(p_1,p) \cdot \xi_{\rm int} \cdot
\frac{V_p^{\rm source} \cdot p_1^2 \Delta p_1 \Delta \xi_1}{V_p^{\rm target} \cdot p^2 \Delta p \Delta \xi}
\]

5. 写入 PETSc 矩阵，最终组装。

### 运行时流程

```
求解器每时间步：
  ├─ Rebuild(t, dt) ── 空操作（矩阵已预计算）
  ├─ BuildJacobian()
  │    └─ SetJacobianBlock(f_re, f_re, jac, x)
  │         └─ SetMatrixElements(jac, rhs)
  │              └─ 对每个径向位置 ir：
  │                   读取 K[ir] 的每一行，
  │                   乘以 scaleFactor × n_tot[ir]，
  │                   用 mat->SetRow() 写入 FVM 矩阵
  │
  └─ 求解线性系统
       └─ SetVectorElements(vec, x)
            └─ 对每个径向位置 ir：
                  Vec = f_re[ir]
                  MatMult(K[ir], f_re_vec, S_vec)
                  vec[ir] += scaleFactor × n_tot[ir] × S_vec
```

#### 关于局部索引

`FVM::Matrix` 的 `SetRow(row, ncols, cols, vals)` 方法会自动加上 `rowOffset` 和 `colOffset`（由求解器的 `SelectSubEquation()` 设置），将块内索引映射到全局矩阵的正确位置。

因此传入的 `row` 和 `cols` 必须使用**块内局部索引**（0 到 NCellsPerRadius-1），而不是全局索引。传入全局索引会导致双重偏移，写入错误的内存位置。

### 雅可比矩阵

`SetJacobianBlock` 提供解析雅可比：

| 求导对象 | 表达式 | 实现 |
|----------|--------|------|
| \(\partial S / \partial f_{\rm re}\) | \(n_{\rm tot} \cdot {\rm scaleFactor} \cdot K\) | 调用 `SetMatrixElements(jac)` |
| \(\partial S / \partial n_{\rm tot}\) | \({\rm scaleFactor} \cdot K \cdot f_{\rm re}\) | 对角矩阵：每径向点一列 |

---

## 注册与接入

### 文件清单

| 文件 | 操作 | 说明 |
|------|------|------|
| `include/DREAM/Equations/Kinetic/AvalancheSourceDirect.hpp` | 新建 | 类声明 |
| `src/Equations/Kinetic/AvalancheSourceDirect.cpp` | 新建 | 完整实现 |
| `src/CMakeLists.txt` | 编辑 | 注册源文件到编译 |
| `include/DREAM/Settings/OptionConstants.enum.hpp` | 编辑 | 添加 `EQTERM_AVALANCHE_MODE_DIRECT = 5` |
| `py/DREAM/Settings/Equations/RunawayElectrons.py` | 编辑 | 添加 `AVALANCHE_MODE_DIRECT = 5`，验证逻辑 |
| `src/Settings/Equations/f_re.cpp` | 编辑 | 在 `ConstructEquation_f_re_kineq()` 中注册 |

### 注册方式

`AvalancheSourceDirect` 通过 `eqn->AddTerm(...)` 直接添加到已有的主算子中：

```cpp
// f_re.cpp: ConstructEquation_f_re_kineq()
FVM::Operator *eqn = ConstructEquation_f_general(...);  // 创建主算子(含电场/碰撞/...)
                                                         // 并在内部注册 eqsys->SetOperator(id_f_re, id_f_re, eqn)

// 在流体源项之后，检查是否需要添加雪崩源
if (ava_mode == EQTERM_AVALANCHE_MODE_DIRECT) {
    eqn->AddTerm(new AvalancheSourceDirect(runawayGrid, unknowns, pCut, -1.0));
}
```

**关键注意**：不能使用 `eqsys->SetOperator(id_f_re, id_f_re, ...)` 注册雪崩源算子，因为 `SetOperator` 内部是 `map[key] = value`，会**覆盖**已有的主算子（包含瞬态项、电场、碰撞、同步辐射等），导致方程系统缺失核心项、雅可比奇异、求解器崩溃。

### Python 接口

```python
import DREAM.Settings.Equations.RunawayElectrons as Runaways

ds.eqsys.n_re.setAvalanche(
    avalanche=Runaways.AVALANCHE_MODE_DIRECT,
    pCutAvalanche=2.0       # 动量下限（必须显式设置）
)
```

### 与其他雪崩模式的关系

| 模式 | 枚举值 | 注册位置 | 算子类型 | 对角/非对角 |
|------|--------|----------|----------|------------|
| NEGLECT | 1 | — | — | — |
| FLUID | 2 | `RunawaySourceTerms.cpp` | 流体生长率 | f_re→n_re 非对角 |
| FLUID_HESSLOW | 3 | `RunawaySourceTerms.cpp` | 流体生长率(Hesslow) | f_re→n_re 非对角 |
| KINETIC (RP) | 4 | `RunawaySourceTerms.cpp` + `f_hot.cpp` | Rosenbluth-Putvinski 解析公式 | f_re→n_re/n_tot/n_i 非对角 |
| **DIRECT** | **5** | **`f_re.cpp`** | **Møller 截面直接数值积分** | **f_re→f_re 对角** |

### 已知问题与注意事项

1. **PETSc 销毁顺序**：`main.cpp` 中 `dream_finalize()`（含 `PetscFinalize()`）在 `delete sim` 之前调用，导致 `AvalancheSourceDirect` 析构函数中的 `MatDestroy` 在 PETSc 已终结后执行。通过 `DeallocateMatrices()` 中检查 `PetscInitialized()` 来保护。
2. **数值保护**：`asin` 参数钳制到 [-1,1]、Vp ≤ 0 跳过、`isfinite` 过滤，防止浮点异常传播到求解器。
3. **计算效率**：PETSc 稀疏矩阵在每次 `SetVectorElements` 和 `SetJacobianBlock` 中都会被遍历，这是一个稠密的动量空间卷积，计算量正比于矩阵非零元数。

### 编译与验证

```bash
cd /data/zhzhou/DREAM/build
cmake .. && make -j$(nproc)
```

验证运行：

```bash
cd /data/zhzhou/DREAM/examples/avalanche_whistler/scripts
./run_DIIID_fre.sh
```
