# Fortran (1.f90) vs CUDA (rk4_test.cu) 数学实现细节对比

## 概述

两个文件均为求解二维 Euler 方程组的 Discontinuous Galerkin (DG) 有限元方法，算例为经典的 Double Mach Reflection (双马赫反射) 问题。Fortran 版本为 P1 阶 MPI 并行实现，CUDA 版本为 P3 阶单 GPU 实现。

### 控制方程

$$
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}(\mathbf{U})}{\partial x} + \frac{\partial \mathbf{G}(\mathbf{U})}{\partial y} = 0
$$

其中守恒变量 $\mathbf{U} = (\rho, \rho u, \rho v, E)^T$，$\gamma = 1.4$。

---

## 1. 多项式阶数与基函数体系

### 1.1 Fortran: P1 阶, 修正正交基

- **多项式阶数**: $k=1$ (P1)
- **基函数总数**: $\dim P_k = \frac{(k+1)(k+2)}{2} = 3$
- **基函数定义**: 采用缩放至参考单元 $[-1,1]^2$ 的修正正交多项式，手动编码

$$
\varphi_1(\xi,\eta) = 1, \quad \varphi_2(\xi,\eta) = \xi, \quad \varphi_3(\xi,\eta) = \eta
$$

- **质量矩阵的对角元** (利用正交性):

$$
m_1 = 1, \quad m_2 = \frac{1}{3}, \quad m_3 = \frac{1}{3}
$$

- **数值积分精度**: 5 点 Gauss-Legendre 求积 (NumGLP=5)，对应 $2 \times 5 - 1 = 9$ 阶代数精度，对 P1 阶体积分远超必要精度
- **面求值点**: 同样为 5 点 Gauss 求积，在每个面上取 5 个点计算数值通量

**关键细节**: Fortran 代码中基函数的导数已经预计算并存储：

$$
\frac{\partial \varphi_1}{\partial \xi} = 0, \quad \frac{\partial \varphi_2}{\partial \xi} = \frac{1}{h_x/2}, \quad \frac{\partial \varphi_3}{\partial \xi} = 0
$$

$$
\frac{\partial \varphi_1}{\partial \eta} = 0, \quad \frac{\partial \varphi_2}{\partial \eta} = 0, \quad \frac{\partial \varphi_3}{\partial \eta} = \frac{1}{h_y/2}
$$

这里 $\frac{1}{h_x/2} = \frac{2}{h_x}$ 是从参考坐标到物理坐标的链式法则因子。

### 1.2 CUDA: P3 阶, Legendre 张量积基

- **多项式阶数**: $k=3$ (P3)
- **基函数总数**: $N_{\text{MODE}} = (k+1)(k+1) = 10$ (采用张量积，非三角形上的缩简基)
- **基函数索引映射**:

$$
\varphi_k(\xi,\eta) = L_{m_k}(\xi) \cdot L_{n_k}(\eta), \quad k = 0,1,\ldots,9
$$

其中 $(m_k, n_k)$ 的映射为:

| k | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| $m_k$ | 0 | 1 | 0 | 2 | 1 | 0 | 3 | 2 | 1 | 0 |
| $n_k$ | 0 | 0 | 1 | 0 | 1 | 2 | 0 | 1 | 2 | 3 |

即按总阶数 $m_k + n_k$ 从小到大排列。

- **Legendre 多项式**:

$$
L_0(x) = 1, \quad L_1(x) = x, \quad L_2(x) = \frac{1}{2}(3x^2 - 1), \quad L_3(x) = \frac{1}{2}(5x^3 - 3x)
$$

- **质量矩阵逆** (利用 Legendre 正交性):

$$
M_{kk}^{-1} = \frac{(2m_k+1)(2n_k+1)}{4}, \quad k = 0,1,\ldots,9
$$

代码实现为 `M_diag_inv_h[k] = 1.0 / ((2.0/(2.0*mk+1.0)) * (2.0/(2.0*nk+1.0)))`。

- **数值积分**: 4 点 Gauss 求积 (N_FACE_PTS=4)，体积分用 $4 \times 4 = 16$ 个张量积点 (N_QUAD=16)

**关键区别**: Fortran 用的是 **三角形上的缩简基** $\varphi_d$（总自由度 $(k+1)(k+2)/2$），而 CUDA 用的是 **矩形上的张量积基**（总自由度 $(k+1)^2$）。这是两种完全不同的基函数体系，导致自由度数量不同：P1 时 Fortran 有 3 个自由度，P3 时 CUDA 有 10 个自由度。

---

## 2. 初始条件的数学表达

### 2.1 激波位置

激波初始位置为过点 $(x_0, y_0) = (1/6, 0)$、法向角 $\theta = \pi/6$ 的直线:

$$
x_s(y) = \frac{1}{6} + \frac{y}{\tan(\pi/3)} = \frac{1}{6} + \frac{y}{\sqrt{3}}
$$

### 2.2 Fortran 的初始条件

```
if (x < 1d0/6d0 + y/3d0**0.5) then  ! 激波后 (post-shock)
    rho = 8,   u = 8.25*cos(pi/6),   v = -8.25*sin(pi/6),   p = 116.5
else                                   ! 激波前 (pre-shock)
    rho = 1.4, u = 0,                 v = 0,                  p = 1
end if
```

激波后速度分量:
- $u_1 = 8.25 \cos(\pi/6) = 8.25 \cdot \frac{\sqrt{3}}{2} \approx 7.1446$
- $v_1 = -8.25 \sin(\pi/6) = -8.25 \cdot \frac{1}{2} = -4.125$

### 2.3 CUDA 的初始条件

```
u = 4.125 * sqrt(3.0) ≈ 7.1446    (与 Fortran 一致)
v = -4.125                          (与 Fortran 一致)
```

数学上完全一致。但 CUDA 直接用数值结果 $4.125\sqrt{3}$ 替代了 $8.25\cos(\pi/6)$，避免了浮点三角函数调用。

### 2.4 初始投影

**Fortran** 采用 $L^2$ 投影，在 Gauss 点上精确积分:

$$
\hat{u}_d = \frac{1}{m_d} \sum_{i_1,j_1} \frac{1}{4} w_{i_1} w_{j_1} \, U(x_c^i + h_x/2 \cdot \lambda_{i_1},\; y_c^j + h_y/2 \cdot \lambda_{j_1}) \, \varphi_d(\lambda_{i_1}, \lambda_{j_1})
$$

**CUDA** 仅保留均值模态（第 0 阶系数），高阶系数初始化为零:

$$
\hat{u}_0^{(v)} = Q_v(x_c, y_c), \quad \hat{u}_k^{(v)} = 0 \;\; (k \geq 1)
$$

这是一个重大区别：Fortran 的初始投影是精确的 $L^2$ 投影，而 CUDA 只用了逐单元的常数近似。这意味着 CUDA 的初始解只有零阶精度，需要若干时间步让高阶模态通过 DG 演化自然生成。

---

## 3. 空间离散：弱形式

### 3.1 DG 弱形式

在单元 $K_{ij}$ 上，试探函数 $\varphi_d$，DG 弱形式为:

$$
\int_{K_{ij}} \frac{\partial \mathbf{U}_h}{\partial t} \varphi_d \, d\mathbf{x} = \underbrace{\int_{K_{ij}} (\mathbf{F} \cdot \nabla_\xi \varphi_d + \mathbf{G} \cdot \nabla_\eta \varphi_d) \, d\mathbf{x}}_{\text{体积分}} - \underbrace{\oint_{\partial K_{ij}} \hat{\mathbf{F}}_n \, \varphi_d \, ds}_{\text{面积分}}
$$

### 3.2 Fortran 的体积分实现

```fortran
! 体积积分：只在 d > 1 时计算（d=1 对应常数基函数，其梯度为零）
du(i,j,k,d,n) = du(i,j,k,d,n) + 0.25*weight(i1)*weight(j1)*
    (Fx(i1,j1,k,n)*phixG(i1,j1,d) + Fy(i1,j1,k,n)*phiyG(i1,j1,d))
```

数学表达:

$$
\text{Vol}_d = \frac{h_x h_y}{4} \sum_{i_1,j_1} w_{i_1} w_{j_1} \left[ F(\mathbf{U}_h(\xi_{i_1},\eta_{j_1})) \frac{\partial \varphi_d}{\partial x} + G(\mathbf{U}_h(\xi_{i_1},\eta_{j_1})) \frac{\partial \varphi_d}{\partial y} \right]
$$

其中 $\frac{1}{4}$ 因子来自面积变换 $dxdy = \frac{h_x}{2}\frac{h_y}{2} d\xi d\eta$。

**注意**: Fortran 代码中 `d > 1` 的判断是因为 $\varphi_1 \equiv 1$ 的梯度为零，因此体积分只对 $d = 2, 3$ 计算。

### 3.3 CUDA 的体积分实现

```c
Vol_Int[v][k] += d_w_quad[q] * (F_val[v] * d_dphi_dr_vol[k][q] * (dy/2.0) +
                                  G_val[v] * d_dphi_ds_vol[k][q] * (dx/2.0));
```

数学表达:

$$
\text{Vol}_{v,k} = \sum_q w_q \left[ F_v(\mathbf{U}_h) \frac{\partial \varphi_k}{\partial \xi} \frac{h_y}{2} + G_v(\mathbf{U}_h) \frac{\partial \varphi_k}{\partial \eta} \frac{h_x}{2} \right]
$$

**关键数学区别**: Fortran 中导数 $\text{phixG} = \frac{\partial \varphi}{\partial x}$ 已经包含了 Jacobian 因子 $\frac{2}{h_x}$，所以求积权重只乘 $\frac{1}{4}$。而 CUDA 中 `d_dphi_dr_vol` 是参考域上的导数 $\frac{\partial \varphi}{\partial \xi}$，不含 Jacobian，因此 Jacobian 因子 $\frac{h_y}{2}$ 和 $\frac{h_x}{2}$ 分别乘在各项上。两种写法数学上等价。

### 3.4 Fortran 的面积分实现

```fortran
! x 方向面积分
du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hx)*weight(j1)*
    (Fxhat(i,j,k,j1,n)*phiGR(j1,d) - Fxhat(i-1,j,k,j1,n)*phiGL(j1,d))

! y 方向面积分
du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hy)*weight(i1)*
    (Fyhat(i,j,k,i1,n)*phiGU(i1,d) - Fyhat(i,j-1,k,i1,n)*phiGD(i1,d))
```

x 方向面积分的数学含义:

$$
\text{Surf}_{x,d} = \frac{1}{h_x/2} \sum_{j_1} w_{j_1} \left[ \hat{F}_x(\xi=+1, \eta_{j_1}) \varphi_d(+1, \eta_{j_1}) - \hat{F}_x(\xi=-1, \eta_{j_1}) \varphi_d(-1, \eta_{j_1}) \right] \cdot \frac{h_y}{2}
$$

其中 $\frac{0.5}{h_x}$ 来源于 $\frac{\partial \varphi}{\partial x}$ 与面法向的缩放（实际上是分部积分后 $1/h_x$ 的因子，乘以参考单元长度 1/2）。

`phiGR(j1,d)` = $\varphi_d(+1, \eta_{j_1})$ 是右面基函数值，`phiGL(j1,d)` = $\varphi_d(-1, \eta_{j_1})$ 是左面基函数值。

### 3.5 CUDA 的面积分实现

```c
Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_B[k][p] * (dx/2.0);  // 底面
Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_T[k][p] * (dx/2.0);  // 顶面
Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_L[k][p] * (dy/2.0);  // 左面
Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_R[k][p] * (dy/2.0);  // 右面
```

注意: 底面法向量 $\mathbf{n} = (0,-1)$，所以 `llf_flux_vector(U_minus, U_p, 0, -1, num_f)` 的第四个参数为 $-1$。

最终 RHS 组装:

$$
\text{RHS}_{v,k} = \frac{M_{kk}^{-1}}{J} (\text{Vol}_{v,k} - \text{Surf}_{v,k}), \quad J = \frac{h_x h_y}{4}
$$

而 Fortran 最后一步是:

$$
du_d = \frac{du_d}{m_d}
$$

即除以质量矩阵对角元 $m_d$。两者数学等价: $\frac{1}{m_d \cdot J}$ vs $\frac{M_{kk}^{-1}}{J}$。

### 3.6 最终除以质量矩阵

**Fortran**:

```fortran
do d = 1,dimPk1
    du(:,:,:,d,:) = du(:,:,:,d,:)/mm(d)
end do
```

**CUDA**:

```c
d_RHS[idx][v][k] = (Vol_Int[v][k] - Surf_Int[v][k]) * d_M_diag_inv[k] / J;
```

两者都实现了 $\hat{u}_d' = M_{dd}^{-1} \times (\text{Vol} - \text{Surf})$，但 CUDA 把 Jacobian $J$ 也一并处理了。

---

## 4. 数值通量

### 4.1 Fortran: LF 与 HLL 双选项

**Lax-Friedrichs (LF) 通量** (`flux_type = 1`):

$$
\hat{\mathbf{F}}_{\text{LF}} = \frac{1}{2}\left[\mathbf{F}(\mathbf{U}_L) + \mathbf{F}(\mathbf{U}_R) - \alpha_{\max}(\mathbf{U}_R - \mathbf{U}_L)\right]
$$

其中 $\alpha_{\max} = \max(|S_R|, |S_L|)$，$S_R = u_n + c$，$S_L = u_n - c$。

**HLL 通量** (`flux_type = 2`):

$$
\hat{\mathbf{F}}_{\text{HLL}} = \begin{cases}
\mathbf{F}_L & \text{if } S_R < 0 \\
\frac{S_R \mathbf{F}_L - S_L \mathbf{F}_R + S_L S_R (\mathbf{U}_R - \mathbf{U}_L)}{S_R - S_L} & \text{if } S_L \leq 0 \leq S_R \\
\mathbf{F}_R & \text{if } S_L > 0
\end{cases}
$$

**波速估计**: 使用两侧状态的最大/最小特征速度:

$$
S_R = \max(u_{n,L} + c_L,\; u_{n,R} + c_R), \quad S_L = \min(u_{n,L} - c_L,\; u_{n,R} - c_R)
$$

### 4.2 CUDA: 仅 LLF

```c
double alpha = fmax(fabs(unL) + cL, fabs(unR) + cR);
flux_res[v] = 0.5 * (flux_n_L + flux_n_R - alpha * (UR[v] - UL[v]));
```

这是 **Local Lax-Friedrichs (LLF/Rusanov)** 通量，与 Fortran 的 LF 基本相同。但 CUDA 的实现是通量形式 $\mathbf{F}\cdot\mathbf{n}$，而 Fortran 分别处理 x 和 y 方向。

**重要区别**: Fortran 的 LF 通量使用 `max(abs(SR), abs(SL))` 作为耗散系数，即取所有特征波速绝对值的上界。CUDA 用的是 `max(|u_n| + c)` 的两侧最大值，这与 `max(|SR|, |SL|)` 在数学上等价（因为 $|u_n - c| \leq |u_n| + c$）。

### 4.3 通量计算中的 Gauss 点数差异

- **Fortran**: 在每条边上取 **5 个 Gauss 点**，在每个点独立计算数值通量
- **CUDA**: 在每条边上取 **4 个 Gauss 点**，同样独立计算

这意味着 Fortran 在通量计算中的分辨率更高，但对于 P1 阶只需要 1 个点就足够精确（边上线性多项式），5 点显然是多余的。CUDA 的 4 点对 P3 阶边上的 4 次多项式（P3 基函数在边上是 3 次，通量是 7 次）仍然是低阶积分。

---

## 5. 时间积分: 10 步四阶强稳定性保持 Runge-Kutta

两个实现采用了相同的非标准 10 步 SSP-RK4 方法，但实现方式不同。

### 5.1 数学形式

设 $\mathbf{L}_h(\mathbf{U})$ 为空间离散算子，10 步方法为:

**阶段 I** (步 1-5): 以 $\Delta t/6$ 为子步长的前向 Euler

$$
\mathbf{U}^{(i)} = \mathbf{U}^{(i-1)} + \frac{\Delta t}{6} \mathbf{L}_h(\mathbf{U}^{(i-1)}), \quad i = 1, 2, 3, 4, 5
$$

每步后施加 jump filter 和保正性限制器。

**中间重组**:

$$
\mathbf{U}_{II} = 0.04\,\mathbf{U}^n + 0.36\,\mathbf{U}^{(5)}
$$

$$
\mathbf{U}^{(5)'} = 15\,\mathbf{U}_{II} - 5\,\mathbf{U}^{(5)}
$$

**阶段 II** (步 6-9): 继续 $\Delta t/6$ 子步长

$$
\mathbf{U}^{(i)} = \mathbf{U}^{(i-1)} + \frac{\Delta t}{6} \mathbf{L}_h(\mathbf{U}^{(i-1)}), \quad i = 6, 7, 8, 9
$$

**最终组合**:

$$
\mathbf{U}^{n+1} = \mathbf{U}_{II} + 0.6\,\mathbf{U}^{(9)} + \frac{\Delta t}{10} \mathbf{L}_h(\mathbf{U}^{(9)})
$$

### 5.2 Fortran 实现

```fortran
uI = uh
uII = uh                    ! 保存 U^n

do i = 1,5
    call Lh                 ! 计算 du = Lh(uh)
    uI = uh + (dt/6d0)*du   ! 前向 Euler 子步
    uh = uI
    call apply_jump_filter_limiter  ! 每步施加限制器
end do

uII = 0.04d0*uII + 0.36d0*uI   ! 重组
uI = 15*uII - 5*uI             ! 第二阶段起始
uh = uI
tRK = tRK - 0.5*dt             ! 时间回退半步

do i = 6,9
    call Lh
    uI = uh + (dt/6d0)*du
    uh = uI
    call apply_jump_filter_limiter
end do

call Lh
uh = uII + 0.6d0*uI + (dt/10d0)*du  ! 最终组合
call apply_jump_filter_limiter
```

**关键**: Fortran 在子步时间上推进 `tRK`，阶段 I 从 $t^n$ 推到 $t^n + \frac{5}{6}\Delta t$，然后回退到 $t^n + \frac{1}{2}\Delta t$，阶段 II 从 $t^n + \frac{1}{2}\Delta t$ 推到 $t^n + \frac{9}{6}\Delta t$。

### 5.3 CUDA 实现

CUDA 使用三个独立 kernel 拆分时间推进:

1. `update_rk4_stage_forward`: 前向 Euler 子步 $\mathbf{U} \mathrel{+}= \frac{\Delta t}{6}\mathbf{L}_h$
2. `update_rk4_intermediate`: 中间重组 $\mathbf{U}_{II} = 0.04\mathbf{U}^n + 0.36\mathbf{U}^{(5)}$
3. `update_rk4_final`: 最终组合 $\mathbf{U}^{n+1} = \mathbf{U}_{II} + 0.6\mathbf{U}^{(9)} + \frac{\Delta t}{10}\mathbf{L}_h(\mathbf{U}^{(9)})$

**重要区别**: CUDA 的 `update_rk4_intermediate` 中:

```c
double uII_new = 0.04 * u_n + 0.36 * u_5_filtered;
d_Mesh[id].U[v][k] = 15.0 * uII_new - 5.0 * u_5_filtered;
```

这里 `u_5_filtered` 是经过 Zhang-Shu 限制器后的状态，而不是未滤波的状态。这确保了进入阶段 II 的数据是物理上合法的。Fortran 中也类似——`uI` 在赋值前已经过 `apply_jump_filter_limiter`。

### 5.4 SSP 性质讨论

这个 10 步方法不是标准的 SSP-RK。标准四阶 SSP 方法需要至少 5 个阶段，且 CFL 系数 $\leq 1.508$。此处的 10 步格式通过在每个子步后施加滤波器和限制器来维持稳定性，但破坏了严格的 SSP 结构。使用的 CFL = 0.75 是经验选取的。

---

## 6. Jump Filter (激波捕捉滤波器)

这是两个实现之间**数学差异最大**的部分。

### 6.1 Fortran: 基于单元面值跳变的二阶导数滤波器

#### Step 1: 计算面上的函数值和导数值

在每条边上取 **2 个端点** (顶点)，对每个守恒变量分量计算:

- **零阶矩跳变** (函数值):

$$
\delta_0^{\max} = \max_{s \in \{\rho, \rho u, \rho v, E\}} \left[\beta_x \sum_{\text{面}} \sum_{v=1}^{2} |U_s^{\text{in}}(v) - U_s^{\text{out}}(v)| + \beta_y \sum_{\text{面}} \sum_{v=1}^{2} |U_s^{\text{in}}(v) - U_s^{\text{out}}(v)|\right]
$$

其中 $\beta_x = |u_1| + c$，$\beta_y = |u_2| + c$ 是局部最大波速。

- **一阶矩跳变** (导数值):

$$
\delta_1^{\max} = \max_{s} \left[\beta_x \cdot 2h_x \sum_{\text{面}} \sum_{v=1}^{2} |\nabla U_s^{\text{in}}(v) - \nabla U_s^{\text{out}}(v)| + \beta_y \cdot 2h_y \sum_{\text{面}} \sum_{v=1}^{2} |\nabla U_s^{\text{in}}(v) - \nabla U_s^{\text{out}}(v)|\right]
$$

注意: 对于 P1 阶，面上的导数只需要考虑 $\partial_x$ 和 $\partial_y$ 分量（没有交叉导数），所以 `phiGR_ver_derx(:,2) = (1,1)` 和 `phiGR_ver_dery(:,2) = (0,0)` 表示 $\varphi_2 = \xi$ 只有 $x$ 方向导数。

#### Step 2: 计算 damping 系数

```fortran
scal = 0.5d0
scal = scal / enthayij              ! 局部缩放因子 = 0.5 / H
damping = delta0max + delta1max
damping = scal * hx * damping / hx   ! = scal * damping = 0.5 * damping / H
```

数学化简后:

$$
\sigma = \frac{0.5}{H}(\delta_0^{\max} + \delta_1^{\max})
$$

其中 $H = (E + p)/\rho$ 是比焓，用作无量纲化因子。

#### Step 3: 指数衰减

```fortran
uhmod(i,j,k,2:3,1:4) = exp(-dt*damping)*uh(i,j,k,2:3,1:4)
```

$$
\hat{u}_d^{(v)} \leftarrow e^{-\sigma \Delta t} \hat{u}_d^{(v)}, \quad d = 2,3;\; v = 1,\ldots,4
$$

**只衰减第 2 和第 3 个基函数系数**（即 $\xi$ 和 $\eta$ 方向的线性模态），不触碰第 1 个（均值）。

**数学意义**: 对 P1 阶，只有 2 个高阶模态（$d=2,3$），滤波器对它们统一施加一个指数衰减因子。这意味着解在激波附近被退化为逐片常数。

### 6.2 CUDA: 基于顶点导数跳变的多阶分离滤波器

#### Step 1: 预计算顶点处的高阶导数

对每个单元的 4 个顶点 $(\pm 1, \pm 1)$，计算 0-3 阶共 10 个偏导数:

$$
\text{deriv}[d][v] = \sum_{k=0}^{9} \hat{u}_k^{(v)} D_d \varphi_k(\xi, \eta), \quad d = 0,1,\ldots,9
$$

其中 $D_d$ 代表 10 种偏导数算子:

| d | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|
| 阶数 | $u$ | $u_\xi$ | $u_\eta$ | $u_{\xi\xi}$ | $u_{\xi\eta}$ | $u_{\eta\eta}$ | $u_{\xi\xi\xi}$ | $u_{\xi\xi\eta}$ | $u_{\xi\eta\eta}$ | $u_{\eta\eta\eta}$ |
| $m$ | 0 | 1 | 1 | 2 | 2 | 2 | 3 | 3 | 3 | 3 |

注意 Jacobian 因子:

$$
\frac{\partial}{\partial x} = \frac{2}{h_x}\frac{\partial}{\partial \xi}, \quad \frac{\partial}{\partial y} = \frac{2}{h_y}\frac{\partial}{\partial \eta}
$$

#### Step 2: 在顶点处计算跳变

对于每个顶点 $\alpha = 0,1,2,3$（分别对应 $(-1,-1), (+1,-1), (+1,+1), (-1,+1)$），计算当前单元与邻居在共享顶点处各阶导数的差:

$$
\text{jump}[\alpha][\text{dir}][d][v] = |D_d U^{(\text{curr})}_v(\text{vertex}_\alpha) - D_d U^{(\text{neigh})}_v(\text{vertex}_\alpha)|
$$

其中 $\text{dir} \in \{X, Y\}$ 表示跳变来自 x 方向邻居还是 y 方向邻居。

具体地，4 个顶点的邻居对应关系:

| 顶点 | X 方向邻居来源 | Y 方向邻居来源 |
|------|--------------|--------------|
| 0 (-1,-1) | 左邻居的顶点 1 | 下邻居的顶点 2 |
| 1 (+1,-1) | 右邻居的顶点 0 | 下邻居的顶点 3 |
| 2 (+1,+1) | 右邻居的顶点 3 | 上邻居的顶点 0 |
| 3 (-1,+1) | 左邻居的顶点 2 | 上邻居的顶点 1 |

#### Step 3: 按导数阶数分组累加

```c
int order_map[10] = {0, 1, 1, 2, 2, 2, 3, 3, 3, 3};
for(int d=0; d<N_DERIV; d++) {
    int m = order_map[d];
    sum_X[m][s] += 0.5*(jump[0][0][d][s] + jump[3][0][d][s] + jump[1][0][d][s] + jump[2][0][d][s]);
    sum_Y[m][s] += 0.5*(jump[0][1][d][s] + jump[1][1][d][s] + jump[3][1][d][s] + jump[2][1][d][s]);
}
```

$$
S_X^{(m)}[s] = \frac{1}{2}\sum_{\alpha=0}^{3} \sum_{d:\text{order}(d)=m} \text{jump}[\alpha][X][d][s]
$$

即将 10 种导数按总阶数 $m = 0,1,2,3$ 分为 4 组，在 x 方向和 y 方向分别求和。

#### Step 4: 计算各阶的 damping 系数

```c
if (m == 0) term_s = (beta_x * sum_X[0][s] + beta_y * sum_Y[0][s]);
else if (m == 1) term_s = (2.0 * dx * beta_x * sum_X[1][s] + 2.0 * dy * beta_y * sum_Y[1][s]);
else if (m == 2) term_s = (6.0 * dx*dx * beta_x * sum_X[2][s] + 6.0 * dy*dy * beta_y * sum_Y[2][s]);
else if (m == 3) term_s = (12.0 * dx*dx*dx * beta_x * sum_X[3][s] + 12.0 * dy*dy*dy * beta_y * sum_Y[3][s]);
term_s /= H;
```

数学表达:

$$
\sigma^{(m)}_s = \frac{1}{H}\left[\beta_x \cdot C_m \cdot h_x^m \cdot S_X^{(m)}[s] + \beta_y \cdot C_m \cdot h_y^m \cdot S_Y^{(m)}[s]\right]
$$

其中系数 $C_m = \{1, 2, 6, 12\}$ 对应 Legendre 多项式的导数归一化常数：

- $C_0 = 1$ (函数值)
- $C_1 = 2$ (一阶导数，因 $L_1'(x) = 1$，参考域长度 2)
- $C_2 = 6$ ($L_2''(x) = 3$，物理域缩放 $4/h_x^2$)
- $C_3 = 12$ ($L_3'''(x) = 15$，物理域缩放 $8/h_x^3$，实际系数 $15 \times \frac{8}{h_x^3} \times h_x^3 = 12$ 近似)

最终:

$$
\sigma^{(m)} = \max_{s} \sigma^{(m)}_s
$$

#### Step 5: 累积 damping 并逐阶施加

```c
d_damp_local[idx][0] = max_contrib[0];
double running_sum = max_contrib[0];
for (int l = 0; l <= 3; l++) {
    running_sum += max_contrib[l];
    d_damp_local[idx][l] = running_sum;
}
```

$$
\sigma_l = \sum_{m=0}^{l} \sigma^{(m)}, \quad l = 0,1,2,3
$$

然后在滤波器 kernel 中:

```c
int l = d_mk_map[k] + d_nk_map[k];  // 基函数的总阶数
double sigma = d_damp_local[id][l];
d_Mesh[id].U[v][k] *= exp(-sigma * dt);
```

$$
\hat{u}_k^{(v)} \leftarrow e^{-\sigma_{m_k+n_k} \cdot \Delta t} \hat{u}_k^{(v)}
$$

### 6.3 两种 Jump Filter 的数学对比

| 特征 | Fortran (P1) | CUDA (P3) |
|------|-------------|-----------|
| 跳变取样位置 | 面上 2 个端点 (顶点) | 4 个顶点处 |
| 跳变阶数 | 0 阶 (值) + 1 阶 (一阶导) | 0-3 阶 (值到三阶导) |
| 跳变量维度 | 2 个点 $\times$ 2 种导数 $\times$ 4 条面 | 4 个顶点 $\times$ 10 种导数 $\times$ 4 条面 |
| 无量纲化因子 | $\text{scal}/H = 0.5/H$ | $1/H$ |
| 衰减粒度 | 统一衰减所有高阶模态 (d=2,3) | 按 $(m_k+n_k)$ 阶数分 4 级衰减 |
| 总阶数 $l=1$ 的模态 | $\sigma_{\text{total}}$ | $\sigma_0 + \sigma_1$ |
| 总阶数 $l=2$ 的模态 | (不存在) | $\sigma_0 + \sigma_1 + \sigma_2$ |
| 总阶数 $l=3$ 的模态 | (不存在) | $\sigma_0 + \sigma_1 + \sigma_2 + \sigma_3$ |

**核心差异**: Fortran 的滤波器将 0 阶和 1 阶跳变合并为单一 damping 值，统一施加于所有高阶模态。CUDA 的滤波器将跳变按阶数分组，对不同总阶数的基函数施加不同强度的衰减——低阶模态衰减较少，高阶模态衰减更多。这更精细地保持了低阶精度同时抑制高阶振荡。

---

## 7. 保正性限制器 (Positivity-Preserving Limiter)

### 7.1 数学原理

Zhang-Shu 保正性限制器保证密度和压力始终为正。核心思想是: 如果单元均值 $\bar{\rho} > 0$ 且 $\bar{p} > 0$，则存在缩放因子 $\theta \in [0,1]$ 使得:

$$
\tilde{\mathbf{U}}_h = \bar{\mathbf{U}} + \theta(\mathbf{U}_h - \bar{\mathbf{U}})
$$

在所有测试点上满足 $\rho > 0$ 和 $p > 0$。

### 7.2 Fortran: GLL + Gauss 测试点

**测试点**: 使用 `phiGLL` 数组，维度为 `(NumGLP, NumGLP, dimPk, 2)`。第二个指标 $d=1,2$ 对应两组不同的测试点集——一组用 Gauss-Lobatto 点 ($\lambda_L$)，另一组用 Gauss 点 ($\lambda$)。对于 NumGLP=5，每组有 $5 \times 5 = 25$ 个测试点，共 50 个。

**密度限制**:

```fortran
eta1 = abs((uh(i,j,k,1,1) - epsilon)/(uh(i,j,k,1,1) - uhGLL(i1,j1,1,d)))
```

$$
\theta_1 = \min_{(i_1,j_1,d)} \left|\frac{\bar{\rho} - \varepsilon}{\bar{\rho} - \rho(\xi_{i_1}, \eta_{j_1}, d)}\right|, \quad \varepsilon = 10^{-13}
$$

```fortran
if (eta < 1) then
    uh(i,j,k,2:dimPk,1) = 0.9*eta*uh(i,j,k,2:dimPk,1)
end if
```

$$
\hat{u}_d^{(\rho)} \leftarrow 0.9\,\theta_1 \hat{u}_d^{(\rho)}, \quad d \geq 2
$$

**注意**: 0.9 的额外因子提供了安全裕度。

**压力限制**: 使用二分法求解从均值到测试点的直线上压力刚好为零的临界参数 $t_q$:

```fortran
call calculate_tq(uh(i,j,k,1,:),uhGLL(i1,j1,:,d),tq,gamma)
```

二分法实现:

$$
t_{q} = \arg\max\{t \in [0,1] : p(\bar{\mathbf{U}} + t(\mathbf{U}_q - \bar{\mathbf{U}})) > 0\}
$$

容差为 $10^{-14}$。

最终:

$$
\hat{u}_d^{(v)} \leftarrow \theta_2 \hat{u}_d^{(v)}, \quad d \geq 2, \;\; \forall v
$$

其中 $\theta_2 = \min_q t_q$。

**关键**: Fortran 中密度限制只修改密度的模态系数，而压力限制修改**所有**守恒变量的模态系数。

### 7.3 CUDA: Gauss + Gauss-Lobatto 组合测试点

**测试点**: 使用 32 个点 (N_ZS_PTS=32)，由 $4 \times 4$ Gauss 点 + $4 \times 4$ Gauss-Lobatto 点组合:

- 前 16 个: $\{GLL_x\} \times \{Gauss_y\}$（x 方向用 Gauss-Lobatto 点包含端点）
- 后 16 个: $\{Gauss_x\} \times \{GLL_y\}$（y 方向用 Gauss-Lobatto 点包含端点）

```c
// 第一组: x 用 Gauss-Lobatto, y 用 Gauss
nodes_GL_h = {-1.0, -0.4472, 0.4472, 1.0}  // 包含端点 ±1
nodes_G_h  = {-0.8611, -0.3400, 0.3400, 0.8611}

// 第二组: x 用 Gauss, y 用 Gauss-Lobatto
```

**密度限制**:

```c
theta1 = (Ubar[0] - eps) / (Ubar[0] - rho_min);
if (theta1 < 0.0) theta1 = 0.0;
if (theta1 > 1.0) theta1 = 1.0;
for (int v = 0; v < NUM_VARS; v++)
    for (int k = 1; k < N_MODE; k++)
        cell->U[v][k] *= theta1;
```

$$
\theta_1 = \min\left(1,\; \max\left(0,\; \frac{\bar{\rho} - \varepsilon}{\bar{\rho} - \rho_{\min}}\right)\right)
$$

**关键差异**: CUDA 中密度限制后修改的是**所有守恒变量**的模态系数，而 Fortran 只修改密度自身。这意味着 CUDA 的密度限制更保守（对所有变量做了同样的缩放）。

**压力限制**: 同样使用二分法，但迭代次数固定为 50 次:

```c
for (int iter = 0; iter < 50; iter++) {
    double t_mid = 0.5 * (t_L + t_R);
    if (p_mid < eps) t_R = t_mid; else t_L = t_mid;
}
```

50 次迭代的精度: $\Delta t < 2^{-50} \approx 10^{-16}$，接近机器精度。

### 7.4 保正性限制器对比总结

| 特征 | Fortran (P1) | CUDA (P3) |
|------|-------------|-----------|
| 测试点数 | 50 (25×2 组) | 32 (16+16 组) |
| 测试点类型 | Gauss + Gauss-Lobatto | Gauss + Gauss-Lobatto |
| 密度限制作用范围 | 仅密度模态 | 所有守恒变量模态 |
| 额外安全因子 | 0.9 (密度限制) | 无 |
| 二分法收敛判据 | $|t_b - t_a| > 10^{-14}$ | 固定 50 次迭代 |
| 压力限制后更新测试点 | 否 | 是 (用 $\theta_1$ 更新后重新计算 $U_{\text{test}}$) |

CUDA 版本在密度限制后更新测试点值再进行压力限制，这是正确做法——确保压力限制基于已经过密度缩放后的解。

---

## 8. CFL 条件与时间步长

### 8.1 Fortran

```fortran
dt = CFL / (alphax/hx + alphay/hy)
```

$$
\Delta t = \frac{\text{CFL}}{\frac{\max|u+c|}{h_x} + \frac{\max|v+c|}{h_y}}
$$

其中 `alphax` 和 `alphay` 是所有单元和所有 Gauss 点上的最大特征波速，通过 MPI 全局归约得到。

### 8.2 CUDA

```c
double max_wave = thrust::transform_reduce(..., WaveSpeedFunctor(d_Mesh), 0.0, thrust::maximum<double>());
double dt = CFL / max_wave;
```

$$
\Delta t = \frac{\text{CFL}}{\max_{ij}\left(\frac{|u_{ij}| + c_{ij}}{h_x} + \frac{|v_{ij}| + c_{ij}}{h_y}\right)}
```

**数学上等价**，但 Fortran 只用均值模态 ($d=1$) 的特征波速，CUDA 也只用均值 `U[v][0]`。

### 8.3 CFL 数

两者均使用 CFL = 0.75。对于此 10 步 RK4 方法，0.75 是一个经验安全值。标准 P1 DG-SLDG 方法的线性稳定 CFL 上界约为 $\frac{1}{2k+1} = \frac{1}{3}$ (P1) 或 $\frac{1}{7}$ (P3)，但加上限制器后可以取更大的 CFL。

---

## 9. 边界条件的数学处理

### 9.1 Fortran: MPI 通信 + 修正模板

Fortran 通过在 `set_bc` 中用 MPI 交换边界单元数据来填充 ghost 区域。边界条件类型:

| 类型码 | 含义 | 数学处理 |
|--------|------|----------|
| 1 | 周期性 | 直接交换 |
| 2 | 自由出流 | $\mathbf{U}_{\text{ghost}} = \mathbf{U}_{\text{boundary}}$ |
| 3 | 入流/反射混合 | 根据坐标和时间设定固定值或反射 |
| 5 | 反射壁面 | 偶/奇延拓 |

**反射壁面 (bc=5)** 的数学处理:

对于壁面法向为 $x$ 方向:

```fortran
call evenex_y(uh(Nx1,ii,k,:,1),uh(Nx,ii,k,:,1))  ! 密度: 偶延拓
call oddex_y(uh(Nx1,ii,k,:,2),uh(Nx,ii,k,:,2))    ! x 动量: 奇延拓
call evenex_y(uh(Nx1,ii,k,:,3),uh(Nx,ii,k,:,3))   ! y 动量: 偶延拓
call evenex_y(uh(Nx1,ii,k,:,4),uh(Nx,ii,k,:,4))   ! 能量: 偶延拓
```

其中 `evenex_y` 和 `oddex_y` 的实现:

```fortran
subroutine evenex_y(a,b)
    a(1) = b(1)    ! 常数项: 偶
    a(2) = 0d0     ! ξ 方向系数: 置零
    a(3) = 0d0     ! η 方向系数: 置零
end subroutine

subroutine oddex_y(a,b)
    a(1) = -b(1)   ! 常数项: 奇 (取反)
    a(2) = 0d0     ! ξ 方向系数: 置零
    a(3) = 0d0     ! η 方向系数: 置零
end subroutine
```

**数学含义**: 对于 P1 阶，ghost 单元的均值等于内部单元的均值（偶）或取反（奇），而所有高阶系数被置零。这等价于将 ghost 单元设为逐片常数。

### 9.2 CUDA: 显式 Ghost Cell 数组

```c
// 反射壁面 (底部, x >= 1/6 部分):
for (int k = 0; k < N_MODE; k++) {
    int nk = d_nk_map[k];  // y 方向 Legendre 阶数
    double sign = (nk % 2 == 0) ? 1.0 : -1.0;
    ghost.U[0][k] =  sign * mesh.U[0][k];  // 密度: 偶
    ghost.U[1][k] =  sign * mesh.U[1][k];  // x 动量: 偶
    ghost.U[2][k] = -sign * mesh.U[2][k];  // y 动量: 奇
    ghost.U[3][k] =  sign * mesh.U[3][k];  // 能量: 偶
}
```

**数学含义**: 对于张量积基 $\varphi_k = L_{m_k}(\xi) L_{n_k}(\eta)$，反射在 $y$ 方向 ($\eta \to -\eta$) 时:

- **偶延拓**: $L_{n_k}(-\eta) = (-1)^{n_k} L_{n_k}(\eta)$，偶数阶不变，奇数阶取反
- **奇延拓**: 同上但整体再取反

因此 ghost 单元系数:

$$
\hat{u}_k^{\text{ghost}} = (-1)^{n_k} \hat{u}_k^{\text{inner}} \quad (\text{偶延拓})
$$

$$
\hat{u}_k^{\text{ghost}} = -(-1)^{n_k} \hat{u}_k^{\text{inner}} \quad (\text{奇延拓})
$$

**关键区别**: Fortran 将 ghost 单元的高阶系数全部置零（退化为常数），而 CUDA 保留了 P3 的完整模态结构。CUDA 的处理在数学上更精确，因为它正确利用了 Legendre 多项式的奇偶性来构造 ghost 单元，而不是简单地截断为逐片常数。

### 9.3 上边界动态激波位置

**Fortran**:

```fortran
if (Xc(ii) < 1d0/6d0 + (1 + 20*tRK)/3d0**0.5) then  ! 激波后
```

激波位置: $x_s(t) = \frac{1}{6} + \frac{1 + 20t}{\sqrt{3}}$

**CUDA**:

```c
double xs_top = (1.0/6.0) + (10.0*current_time + 0.5) / sin60;
```

激波位置: $x_s(t) = \frac{1}{6} + \frac{10t + 0.5}{\sin(\pi/3)} = \frac{1}{6} + \frac{10t + 0.5}{\sqrt{3}/2}$

**差异**: Fortran 的激波速度参数化为 $20/\sqrt{3}$，CUDA 为 $10/(\sqrt{3}/2) = 20/\sqrt{3}$，数学上一致。但 Fortran 的偏移量为 $1/\sqrt{3}$，CUDA 为 $0.5/(\sqrt{3}/2) = 1/\sqrt{3}$，也一致。

### 9.4 下边界条件

**Fortran** (bcD=3):

```fortran
if (Xc(ii) < 1d0/6d0) then
    uh(ii,0,:,:,:) = uh(ii,1,:,:,:)    ! x < 1/6: 自由出流
else
    call evenex_y(...)                   ! x >= 1/6: 反射壁面
end if
```

**CUDA**:

```c
if (xc_bottom < 1.0/6.0) {
    ghost = Q1;  // 固定入流状态
} else {
    // 反射壁面
}
```

**差异**: Fortran 在 $x < 1/6$ 部分用**自由出流**（复制内部单元值），CUDA 用**固定入流**（设为激波后状态）。数学上，CUDA 的处理更符合物理——下边界 $x < 1/6$ 部分应该始终保持激波后状态。

---

## 10. 计算域与网格

| 参数 | Fortran | CUDA |
|------|---------|------|
| x 范围 | $[0, 3.2]$ | $[0, 4.0]$ |
| y 范围 | $[0, 1.0]$ | $[0, 1.0]$ |
| 全局网格 | $768 \times 240$ | $960 \times 240$ |
| 单元尺寸 $h_x$ | $3.2/768 \approx 0.004167$ | $4.0/960 \approx 0.004167$ |
| 单元尺寸 $h_y$ | $1.0/240 \approx 0.004167$ | $1.0/240 \approx 0.004167$ |
| 有效自由度/变量 | $768 \times 240 \times 3 = 552{,}960$ | $960 \times 240 \times 10 = 2{,}304{,}000$ |

虽然网格尺寸 $h_x, h_y$ 近似相同，但 CUDA 版本的 P3 阶意味着每个单元有 10 个自由度（vs Fortran 的 3 个），有效自由度约为 Fortran 的 4.2 倍，对应更高的空间分辨率。

---

## 11. 特征结构 (Riemann 求解器的辅助计算)

### 11.1 Fortran: 完整特征分解

Fortran 实现了 `compute_Rinv` 子程序，计算 2D Euler 方程的完整右特征向量矩阵 $\mathbf{R}$ 及其逆 $\mathbf{R}^{-1}$:

$$
\mathbf{R} = \begin{pmatrix} 1 & 1 & 1 & 0 \\ u - c n_x & u & u + c n_x & n_y \\ v - c n_y & v & v + c n_y & -n_x \\ H - c u_n & \frac{1}{2}|\mathbf{u}|^2 & H + c u_n & u_y n_x - u_x n_y \end{pmatrix}
$$

其中 $u_n = \mathbf{u} \cdot \mathbf{n}$，$H = (E+p)/\rho$。这在当前代码中虽然定义了但未被 LF/HLL 通量直接使用（这些通量不需要完整的特征分解），可能用于其他限制器或未来扩展。

### 11.2 CUDA: 无特征分解

CUDA 只需要计算最大波速 $\alpha = |u_n| + c$，不需要完整的特征向量。

---

## 12. 总结: 核心数学差异

1. **基函数体系完全不同**: Fortran 用三角域缩简基（3 个 P1 模态），CUDA 用矩形张量积 Legendre 基（10 个 P3 模态）。

2. **初始投影精度不同**: Fortran 做 $L^2$ 投影（精确），CUDA 只用逐单元常数近似。

3. **Jump Filter 的精细程度不同**: Fortran 合并 0-1 阶跳变为单一 damping，统一施加；CUDA 将 0-3 阶跳变分组，按基函数总阶数分 4 级施加不同强度衰减。

4. **保正性限制器的作用域不同**: Fortran 的密度限制只修改密度模态（带 0.9 安全因子），CUDA 修改所有守恒变量模态（无安全因子）。

5. **反射边界条件的精度不同**: Fortran 将 ghost 单元退化为逐片常数，CUDA 保留完整高阶模态结构并利用 Legendre 奇偶性。

6. **数值通量选项不同**: Fortran 支持 LF 和 HLL 两种通量，CUDA 仅使用 LLF。

7. **网格计算域不同**: Fortran 为 $[0,3.2]\times[0,1]$，CUDA 为 $[0,4.0]\times[0,1]$，但单元尺寸相近。
