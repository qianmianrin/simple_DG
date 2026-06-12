#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>  // Windows 下需要这个头文件
#define M_PI 3.14159265358979323846
/* ============================================================
 * 基本全局参数
 * ============================================================ */
#define Nx 50
#define Ny 50
#define xL 0.0
#define xR 1.0
#define yL 0
#define yR 1.0
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma 1.4
#define NUM_VARS 4 

/* 时间推进参数 */
#define CFL   0.1
#define T_END 0.7

/* ============================================================
 * 数据结构: U 中存储的是 Modal 系数 (阶数 k=0...8)
 * ============================================================ */
typedef struct {
    double U[NUM_VARS][9]; 
    double face_top[NUM_VARS][3];
    double face_bottom[NUM_VARS][3];
    double face_left[NUM_VARS][3];
    double face_right[NUM_VARS][3];
} Element;

//对以后的WENO limiter，也许需要不止一层Ghost Cell？
Element Mesh[Nx * Ny];
Element Ghost_bottom[Nx]; 
Element Ghost_top[Nx];    
Element Ghost_left[Ny];   
Element Ghost_right[Ny];  

/* ============================================================
 * 3 点 Gauss-Legendre 节点与权重 (用于积分非线性通量)
 * ============================================================ */
static const double nodes1D[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

double r_quad[9], s_quad[9], w_quad[9]; 

/* ============================================================
 * 预计算矩阵 (Modal DG 特化)
 * ============================================================ */
#define N9 9
double M_diag_inv[N9];       // 对角质量矩阵之逆
double phi_vol[N9][N9];      // [mode_k][quad_q] : 体积分高斯点处的基函数值
double dphi_dr_vol[N9][N9];  // [mode_k][quad_q] : 对 xi 的导数
double dphi_ds_vol[N9][N9];  // [mode_k][quad_q] : 对 eta 的导数
double phi_face_T[N9][3];    // [mode_k][face_quad_p] : 顶面高斯点处的基函数值
double phi_face_B[N9][3];
double phi_face_L[N9][3];
double phi_face_R[N9][3];

/* ============================================================
 * Legendre 基函数定义
 * ============================================================ */
static inline double legendre_1d(int i, double x) {
    if (i == 0) return 1.0;
    if (i == 1) return x;
    if (i == 2) return 0.5 * (3.0 * x * x - 1.0);
    return 0.0;
}

static inline double d_legendre_1d(int i, double x) {
    if (i == 0) return 0.0;
    if (i == 1) return 1.0;
    if (i == 2) return 3.0 * x;
    return 0.0;
}

void init_quadrature(void) {
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) {
            int idx     = j * 3 + i;
            r_quad[idx] = nodes1D[i];
            s_quad[idx] = nodes1D[j];
            w_quad[idx] = weights1D[i] * weights1D[j];
        }
}

void precompute_matrices(void) {
    for (int k = 0; k < N9; k++) {
        int mk = k % 3; // xi 方向多项式阶数
        int nk = k / 3; // eta 方向多项式阶数
        
        // 正交质量矩阵对角元素: M_kk = \int \int P_mk^2 P_nk^2 d\xi d\eta
        double m_val = (2.0 / (2.0 * mk + 1.0)) * (2.0 / (2.0 * nk + 1.0));
        M_diag_inv[k] = 1.0 / m_val;

        // 体积分配置点
        for (int q = 0; q < N9; q++) {
            double r = r_quad[q], s = s_quad[q];
            phi_vol[k][q]     = legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_dr_vol[k][q] = d_legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_ds_vol[k][q] = legendre_1d(mk, r) * d_legendre_1d(nk, s);
        }

        // 面积分配置点
        for (int p = 0; p < 3; p++) {
            double np = nodes1D[p];
            phi_face_T[k][p] = legendre_1d(mk, np) * legendre_1d(nk, 1.0);
            phi_face_B[k][p] = legendre_1d(mk, np) * legendre_1d(nk, -1.0);
            phi_face_L[k][p] = legendre_1d(mk, -1.0) * legendre_1d(nk, np);
            phi_face_R[k][p] = legendre_1d(mk, 1.0) * legendre_1d(nk, np);
        }
    }
}

/* ============================================================
 * 欧拉方程物理通量与状态方程
 * ============================================================ */
static inline double calc_pressure(double rho, double rhou, double rhov, double E) {
    return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}

static inline void euler_flux(double U[NUM_VARS], double F[NUM_VARS], double G[NUM_VARS]) {
    double rho = U[0], rhou = U[1], rhov = U[2], E = U[3];
    double u = rhou / rho, v = rhov / rho;
    double p = calc_pressure(rho, rhou, rhov, E);

    F[0] = rhou;       F[1] = rhou * u + p; F[2] = rhou * v;       F[3] = u * (E + p);
    G[0] = rhov;       G[1] = rhou * v;     G[2] = rhov * v + p;   G[3] = v * (E + p);
}

static inline double max_speed() { 

   
    double max_v = 0.0;
    for (int i = 0; i < Nx * Ny; i++) {
        // 利用第 0 个模态 (单元均值) 来估算波速
        double rho = Mesh[i].U[0][0];
        double rhou = Mesh[i].U[1][0];
        double rhov = Mesh[i].U[2][0];
        double E = Mesh[i].U[3][0];
        
        double u = rhou / rho;
        double v = rhov / rho;
        double p = calc_pressure(rho, rhou, rhov, E);
        double c = sqrt(gamma * p / rho);
        
        max_v = fmax(max_v, fabs(u) + c);
        max_v = fmax(max_v, fabs(v) + c);
    }
    return (max_v < 1e-9) ? 1.0 : max_v; 
}

static inline void lf_flux_vector(double UL[NUM_VARS], double UR[NUM_VARS], double nx, double ny, double flux_res[NUM_VARS]) {
    double FL[NUM_VARS], GL[NUM_VARS], FR[NUM_VARS], GR[NUM_VARS];
    euler_flux(UL, FL, GL); euler_flux(UR, FR, GR);

    double pL = calc_pressure(UL[0], UL[1], UL[2], UL[3]);
    double cL = sqrt(gamma * pL / UL[0]);
    double alphaL = fabs((UL[1]/UL[0])*nx + (UL[2]/UL[0])*ny) + cL;

    double pR = calc_pressure(UR[0], UR[1], UR[2], UR[3]);
    double cR = sqrt(gamma * pR / UR[0]);
    double alphaR = fabs((UR[1]/UR[0])*nx + (UR[2]/UR[0])*ny) + cR;

    double alpha = fmax(alphaL, alphaR);

    for(int v = 0; v < NUM_VARS; v++) 
        flux_res[v] = 0.5 * (FL[v]*nx + GL[v]*ny + FR[v]*nx + GR[v]*ny) - 0.5 * alpha * (UR[v] - UL[v]);
}

/* ============================================================
 * 初始条件: L2 投影到 Modal 空间
 * ============================================================ */
void init_condition(void) {
    double Lx = xR - xL;
    double Ly = yR - yL;
    const double rho_0 = 1.0;
    const double epsilon = 0.2;
    const double u_0 = 1.0;
    const double v_0 = 0.5;
    const double p_0 = 1.0;

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];
            double xc = xL + (ii + 0.5) * dx;
            double yc = yL + (jj + 0.5) * dy;

            // 初始化所有模态系数为 0
            for(int v = 0; v < NUM_VARS; v++)
                for(int k = 0; k < N9; k++) cell->U[v][k] = 0.0;

            // 对每个模态 k 进行 L2 投影积分
            for (int k = 0; k < N9; k++) {
                // 存储 4 个守恒变量在当前模态下的积分累加值
                double sum_Qphi[NUM_VARS] = {0.0};

                for (int q = 0; q < N9; q++) {
                    // 1. 计算高斯点的物理坐标
                    double x_phys = xc + (dx / 2.0) * r_quad[q];
                    double y_phys = yc + (dy / 2.0) * s_quad[q];

                    // 2. 计算该点的初始物理量 (与你的解析解保持一致)
                    double rho = rho_0 + epsilon * sin(2.0 * M_PI * x_phys / Lx) * sin(2.0 * M_PI * y_phys / Ly);
                    double u = u_0;
                    double v = v_0;
                    double p = p_0;

                    // 3. 转换为守恒变量
                    double Q[NUM_VARS];
                    Q[0] = rho;
                    Q[1] = rho * u;
                    Q[2] = rho * v;
                    Q[3] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);

                    // 4. 累加积分: Q * phi_k * weight
                    // 注意：这里是在标准单元 [-1,1] 上积分，Jacobian 已被 M_diag_inv 抵消
                    for (int v = 0; v < NUM_VARS; v++) {
                        sum_Qphi[v] += w_quad[q] * Q[v] * phi_vol[k][q];
                    }
                }

                // 5. 乘以质量矩阵之逆，得到模态系数
                // 这里除以 4.0 是因为标准单元面积是 4，而 M_diag 计算的是单位化的基函数积分
                for (int v = 0; v < NUM_VARS; v++) {
                    cell->U[v][k] = sum_Qphi[v] * M_diag_inv[k] * (1.0 / 1.0); 
                    // 注：如果你 precompute_matrices 里的 M_diag 已经是 \int\int phi^2 dxi deta
                    // 那么这里直接乘 M_diag_inv[k] 即可。
                }
            }
        }
    }
}

/* ============================================================
 * 边界值插值（Modal展开 -> 面 Gauss 点物理值）
 * ============================================================ */
void compute_face_values(Element *cell) {
    for (int var = 0; var < NUM_VARS; var++) {
        for (int p = 0; p < 3; p++) {
            double vt = 0.0, vb = 0.0, vl = 0.0, vr = 0.0;
            for (int k = 0; k < N9; k++) {
                double mode_val = cell->U[var][k];
                vt += phi_face_T[k][p] * mode_val;
                vb += phi_face_B[k][p] * mode_val;
                vl += phi_face_L[k][p] * mode_val;
                vr += phi_face_R[k][p] * mode_val;
            }
            cell->face_top[var][p]    = vt;
            cell->face_bottom[var][p] = vb;
            cell->face_left[var][p]   = vl;
            cell->face_right[var][p]  = vr;
        }
    }
}

void boundary_value(void) {
    for (int idx = 0; idx < Nx * Ny; idx++) compute_face_values(&Mesh[idx]);
}

/* Ghost Cell  */

/* ============================================================
 * 修改后的边界条件处理 (基于面上 Gauss 点坐标判定)
 * ============================================================ */

void apply_ghost_cells(void) {
    // 垂直方向周期性 (Bottom <-> Top)
    for (int i = 0; i < Nx; i++) {
        for (int var = 0; var < NUM_VARS; var++) {
            for (int k = 0; k < 3; k++) {
                // 底部的 Ghost 得到顶部的面值
                Ghost_bottom[i].face_top[var][k] = Mesh[(Ny - 1) * Nx + i].face_top[var][k];
                // 顶部的 Ghost 得到底部的面值
                Ghost_top[i].face_bottom[var][k] = Mesh[0 * Nx + i].face_bottom[var][k];
            }
        }
    }
    // 水平方向周期性 (Left <-> Right)
    for (int j = 0; j < Ny; j++) {
        for (int var = 0; var < NUM_VARS; var++) {
            for (int k = 0; k < 3; k++) {
                Ghost_left[j].face_right[var][k] = Mesh[j * Nx + (Nx - 1)].face_right[var][k];
                Ghost_right[j].face_left[var][k] = Mesh[j * Nx + 0].face_left[var][k];
            }
        }
    }
}

static void get_neighbor_face(int ii, int jj, int face, double u_plus[NUM_VARS][3]) {
    for(int var = 0; var < NUM_VARS; var++) {
        if (face == 0) { 
            if (jj == 0) for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_bottom[ii].face_top[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[(jj - 1) * Nx + ii].face_top[var][k];
        } else if (face == 1) { 
            if (jj == Ny - 1) for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_top[ii].face_bottom[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[(jj + 1) * Nx + ii].face_bottom[var][k];
        } else if (face == 2) { 
            if (ii == 0) for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_left[jj].face_right[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[jj * Nx + (ii - 1)].face_right[var][k];
        } else { 
            if (ii == Nx - 1) for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_right[jj].face_left[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[jj * Nx + (ii + 1)].face_left[var][k];
        }
    }
}

/* ============================================================
 * 限制器 (Modal TVB 极简化)
 * ============================================================ */
static double minmod_tvb(double a, double b, double c, double h) {
    static const double M_TVB = 0.0;
    if (fabs(a) <= M_TVB * h * h) return a;
    if (a > 0.0 && b > 0.0 && c > 0.0) return fmin(a, fmin(b, c));
    if (a < 0.0 && b < 0.0 && c < 0.0) return fmax(a, fmax(b, c));
    return 0.0;
}

static void build_eigen_x(double rho, double u, double v, double p, double R[NUM_VARS][NUM_VARS], double L[NUM_VARS][NUM_VARS]) {
    double c = sqrt(gamma * p / rho), H = (p/(gamma-1.0) + 0.5*rho*(u*u+v*v) + p) / rho;
    double q2 = 0.5*(u*u + v*v), c2 = c*c, gm1 = gamma - 1.0, a1 = gm1 / (2.0*c2), b = gm1 / c2;
    R[0][0]=1.0; R[1][0]=u-c; R[2][0]=v; R[3][0]=H-u*c; R[0][1]=1.0; R[1][1]=u; R[2][1]=v; R[3][1]=q2;
    R[0][2]=0.0; R[1][2]=0.0; R[2][2]=1.0; R[3][2]=v; R[0][3]=1.0; R[1][3]=u+c; R[2][3]=v; R[3][3]=H+u*c;
    L[0][0]= a1*q2 + u/(2.0*c); L[0][1]=-a1*u - 1.0/(2.0*c); L[0][2]=-a1*v; L[0][3]=a1;
    L[1][0]= 1.0 - b*q2; L[1][1]= b*u; L[1][2]= b*v; L[1][3]=-b;
    L[2][0]=-v; L[2][1]= 0.0; L[2][2]= 1.0; L[2][3]=0.0;
    L[3][0]= a1*q2 - u/(2.0*c); L[3][1]=-a1*u + 1.0/(2.0*c); L[3][2]=-a1*v; L[3][3]=a1;
}

static void build_eigen_y(double rho, double u, double v, double p, double R[NUM_VARS][NUM_VARS], double L[NUM_VARS][NUM_VARS]) {
    double c = sqrt(gamma * p / rho), H = (p/(gamma-1.0) + 0.5*rho*(u*u+v*v) + p) / rho;
    double q2 = 0.5*(u*u + v*v), c2 = c*c, gm1 = gamma - 1.0, a1 = gm1 / (2.0*c2), b = gm1 / c2;
    R[0][0]=1.0; R[1][0]=u; R[2][0]=v-c; R[3][0]=H-v*c; R[0][1]=1.0; R[1][1]=u; R[2][1]=v; R[3][1]=q2;
    R[0][2]=0.0; R[1][2]=1.0; R[2][2]=0.0; R[3][2]=u; R[0][3]=1.0; R[1][3]=u; R[2][3]=v+c; R[3][3]=H+v*c;
    L[0][0]= a1*q2 + v/(2.0*c); L[0][1]=-a1*u; L[0][2]=-a1*v - 1.0/(2.0*c); L[0][3]=a1;
    L[1][0]= 1.0 - b*q2; L[1][1]= b*u; L[1][2]= b*v; L[1][3]=-b;
    L[2][0]=-u; L[2][1]= 1.0; L[2][2]= 0.0; L[2][3]=0.0;
    L[3][0]= a1*q2 - v/(2.0*c); L[3][1]=-a1*u; L[3][2]=-a1*v + 1.0/(2.0*c); L[3][3]=a1;
}

static inline double calc_pressure_safe(double rho, double rhou, double rhov, double E) {
    if (rho <= 0.0) return -1.0; 
    return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}
void apply_positivity_limiter(void) {
    const double eps = 1e-13;

    for (int idx = 0; idx < Nx * Ny; idx++) {
        Element *cell = &Mesh[idx];

        /* 真正的单元均值（mode 0，绝对不能被 limiter 修改） */
        double Ubar[NUM_VARS];
        for (int v = 0; v < NUM_VARS; v++) Ubar[v] = cell->U[v][0];

        /* 防止均值本身非物理（守恒性保证，正常不应触发） */
        if (Ubar[0] < eps) {
            for (int k = 1; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++) cell->U[v][k] = 0.0;
            continue;
        }

        /* -------- 收集所有测试点的多项式点值 -------- */
        /* 使用体积 Gauss 点(9) + 四面各3个 Gauss 点(12) = 21 点 */
        double U_test[21][NUM_VARS];
        int pt = 0;

        for (int q = 0; q < N9; q++, pt++) {
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_vol[k][q];
        }
        for (int q = 0; q < 3; q++) {
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_T[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_B[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_L[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N9; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_R[k][q];
            pt++;
        }
        /* pt == 21 */

        /* ======== 第一步：密度保正 ======== */
        /*   theta1 = min(1, (rho_bar - eps) / (rho_bar - rho_min))  */
        double rho_min = U_test[0][0];
        for (int i = 1; i < 21; i++)
            if (U_test[i][0] < rho_min) rho_min = U_test[i][0];

        double theta1 = 1.0;
        if (rho_min < eps) {
            theta1 = (Ubar[0] - eps) / (Ubar[0] - rho_min);
            if (theta1 < 0.0) theta1 = 0.0;
            if (theta1 > 1.0) theta1 = 1.0;

            /* 只缩放 k>=1 的高阶模态；mode 0 (均值) 严格不变 */
            for (int v = 0; v < NUM_VARS; v++)
                for (int k = 1; k < N9; k++)
                    cell->U[v][k] *= theta1;

            /* 同步更新 U_test（用于第二步） */
            for (int i = 0; i < 21; i++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[i][v] = Ubar[v] + theta1 * (U_test[i][v] - Ubar[v]);
        }

        /* ======== 第二步：压力保正 ======== */
        /*   对每个测试点 i，寻找 t in [0,1] 使 p(Ubar + t*(U_test[i]-Ubar)) = eps
             theta2 = min over all i of that t                                     */
        double theta2 = 1.0;
        for (int i = 0; i < 21; i++) {
            double p_i = calc_pressure_safe(U_test[i][0], U_test[i][1],
                                            U_test[i][2], U_test[i][3]);
            if (p_i < eps) {
                /* p(Ubar) > eps 由均值保证；二分法找 t */
                double t_L = 0.0, t_R = 1.0;
                for (int iter = 0; iter < 50; iter++) {
                    double t_mid = 0.5 * (t_L + t_R);
                    double Ut[NUM_VARS];
                    for (int v = 0; v < NUM_VARS; v++)
                        Ut[v] = Ubar[v] + t_mid * (U_test[i][v] - Ubar[v]);
                    double p_mid = calc_pressure_safe(Ut[0], Ut[1], Ut[2], Ut[3]);
                    if (p_mid < eps) t_R = t_mid;
                    else             t_L = t_mid;
                }
                double t_star = t_L;   /* 取保守侧（压力>eps的最大t） */
                if (t_star < theta2) theta2 = t_star;
            }
        }

        /* 只缩放 k>=1 的高阶模态；mode 0 严格不变 */
        if (theta2 < 1.0)
            for (int v = 0; v < NUM_VARS; v++)
                for (int k = 1; k < N9; k++)
                    cell->U[v][k] *= theta2;
    }
}
/* ============================================================
 * 弱形式 DG 右端项 (Weak Form)
 * ============================================================ */
void compute_rhs(double RHS[Nx * Ny][NUM_VARS][N9]) {
    boundary_value();
    apply_ghost_cells();

    double J = (dx * dy) / 4.0; 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double u_plus[NUM_VARS][3];

            double Vol_Int[NUM_VARS][N9] = {{0}};
            double Surf_Int[NUM_VARS][N9] = {{0}};

            /* --- 步骤 A: 体积分贡献 --- */
            // \iint (F * dphi_dx + G * dphi_dy) dxdy
            for (int q = 0; q < N9; q++) {
                double U_phys[NUM_VARS] = {0};
                for (int k = 0; k < N9; k++)
                    for(int v=0; v<NUM_VARS; v++) U_phys[v] += cell->U[v][k] * phi_vol[k][q];
                
                double F_val[NUM_VARS], G_val[NUM_VARS];
                euler_flux(U_phys, F_val, G_val);

                for (int k = 0; k < N9; k++) {
                    for(int v=0; v<NUM_VARS; v++) {
                        // dphi/dx = dphi/dr * (2/dx), dV = (dx*dy/4) * drds -> combined is (dy/2)
                        Vol_Int[v][k] += w_quad[q] * ( F_val[v] * dphi_dr_vol[k][q] * (dy / 2.0) + 
                                                       G_val[v] * dphi_ds_vol[k][q] * (dx / 2.0) );
                    }
                }
            }

            /* --- 步骤 B: 边界积分贡献 --- */
            double U_minus[NUM_VARS], U_p[NUM_VARS], num_f[NUM_VARS];

            // 1. Bottom
            get_neighbor_face(ii, jj, 0, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_bottom[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 0, -1, num_f);
                for(int k=0; k<N9; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_B[k][p] * (dx / 2.0);
            }
            // 2. Top
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 0, 1, num_f);
                for(int k=0; k<N9; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_T[k][p] * (dx / 2.0);
            }
            // 3. Left
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, -1, 0, num_f);
                for(int k=0; k<N9; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_L[k][p] * (dy / 2.0);
            }
            // 4. Right
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 1, 0, num_f);
                for(int k=0; k<N9; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_R[k][p] * (dy / 2.0);
            }

            /* --- 步骤 C: 组合右端项 dU/dt = M^-1 (Vol - Surf) / J --- */
            for(int v=0; v<NUM_VARS; v++) {
                for (int k = 0; k < N9; k++) {
                    RHS[idx][v][k] = (Vol_Int[v][k] - Surf_Int[v][k]) * M_diag_inv[k] / J;
                }
            }
        }
    }
}

/* ============================================================
 * 时间推进 (完全不变)
 * ============================================================ */
static double compute_dt(void) { return CFL * fmin(dx, dy) / max_speed(); }

static double U0[Nx * Ny][NUM_VARS][N9], U1[Nx * Ny][NUM_VARS][N9], U2[Nx * Ny][NUM_VARS][N9], RHS_buf[Nx * Ny][NUM_VARS][N9];
void rk3_step(double dt, int nit) {
    int total = Nx * Ny;
    for (int i = 0; i < total; i++) memcpy(U0[i], Mesh[i].U, NUM_VARS * N9 * sizeof(double));

    // --- 第一步 ---
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++) for (int v = 0; v < NUM_VARS; v++) for (int q = 0; q < N9; q++)
        U1[i][v][q] = U0[i][v][q] + dt * RHS_buf[i][v][q];
    
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U1[i], NUM_VARS * N9 * sizeof(double));

    apply_positivity_limiter(); // <--- 嵌入位置
    for (int i = 0; i < total; i++) memcpy(U1[i], Mesh[i].U, NUM_VARS * N9 * sizeof(double));

    // --- 第二步 ---
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++) for (int v = 0; v < NUM_VARS; v++) for (int q = 0; q < N9; q++)
        U2[i][v][q] = 0.75 * U0[i][v][q] + 0.25 * (U1[i][v][q] + dt * RHS_buf[i][v][q]);
    
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U2[i], NUM_VARS * N9 * sizeof(double));

    apply_positivity_limiter(); // <--- 嵌入位置
    for (int i = 0; i < total; i++) memcpy(U2[i], Mesh[i].U, NUM_VARS * N9 * sizeof(double));

    // --- 第三步 ---
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++) for (int v = 0; v < NUM_VARS; v++) for (int q = 0; q < N9; q++)
        Mesh[i].U[v][q] = (1.0 / 3.0) * U0[i][v][q] + (2.0 / 3.0) * (U2[i][v][q] + dt * RHS_buf[i][v][q]);
    

    apply_positivity_limiter(); // <--- 嵌入位置
}

void output_results(double t) {
    FILE *fp = fopen("result.dat", "w");
    if (fp == NULL) return;
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Rho\", \"U\", \"V\", \"P\"\n");

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double xc = (ii + 0.5) * dx, yc = (jj + 0.5) * dy;

            // 为了与 Nodal 保持相同格式输出，需将 Modal 系数投影回原本的 Gauss 点
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double r = nodes1D[j], s = nodes1D[i];
                    double x_phys = xc + (dx / 2.0) * r, y_phys = yc + (dy / 2.0) * s;
                    
                    double U_phys[NUM_VARS] = {0};
                    for (int k = 0; k < N9; k++) {
                        double phi_val = legendre_1d(k % 3, r) * legendre_1d(k / 3, s);
                        for (int v = 0; v < NUM_VARS; v++) U_phys[v] += cell->U[v][k] * phi_val;
                    }

                    double rho = U_phys[0], u = U_phys[1]/rho, v = U_phys[2]/rho;
                    double p = calc_pressure(rho, U_phys[1], U_phys[2], U_phys[3]);
                    fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", x_phys, y_phys, rho, u, v, p);
                }
            }
        }
    }
    fclose(fp);
    printf("Saved results to result.dat\n");
}

/* ============================================================
 * 误差计算: 计算密度 Rho 的 L1 和 L2 范数误差
 * ============================================================ */
void compute_error(double t_final) {
    double error_L1 = 0.0;
    double error_L2 = 0.0;
    double total_area = 0.0;

    // 初始条件中的参数 (须与 init_condition 保持一致)
    const double rho_0 = 1.0;
    const double epsilon = 0.2;
    const double u_0 = 1.0;
    const double v_0 = 0.5;
    const double Lx = xR - xL;
    const double Ly = yR - yL;

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            
            // 单元中心
            double xc = xL + (ii + 0.5) * dx;
            double yc = yL + (jj + 0.5) * dy;

            // 在单元内使用 3x3 高斯积分点计算误差积分
            for (int q = 0; q < N9; q++) {
                double r = r_quad[q];
                double s = s_quad[q];
                double w = w_quad[q];

                // 高斯点的物理坐标
                double x_phys = xc + (dx / 2.0) * r;
                double y_phys = yc + (dy / 2.0) * s;

                // 1. 获取数值解 rho_num (从 Modal 系数投影)
                double rho_num = 0.0;
                for (int k = 0; k < N9; k++) {
                    rho_num += cell->U[0][k] * phi_vol[k][q];
                }

                // 2. 计算解析解 rho_exact (考虑周期性平移)
                // 使用 fmod 保证坐标回到 [xL, xR] 范围内
                double x_shift = fmod(x_phys - u_0 * t_final - xL, Lx);
                if (x_shift < 0) x_shift += Lx;
                x_shift += xL;

                double y_shift = fmod(y_phys - v_0 * t_final - yL, Ly);
                if (y_shift < 0) y_shift += Ly;
                y_shift += yL;

                double rho_exact = rho_0 + epsilon * sin(2.0 * M_PI * x_shift / Lx) * sin(2.0 * M_PI * y_shift / Ly);

                // 3. 累加误差
                double diff = fabs(rho_num - rho_exact);
                double dA = w * (dx * dy / 4.0); // 雅可比行列式下的积分权重

                error_L1 += diff * dA;
                error_L2 += diff * diff * dA;
                total_area += dA;
            }
        }
    }

    error_L2 = sqrt(error_L2);

    printf("\n================ error analysis (Time: %.2f) ================\n", t_final);
    printf("Mesh Grid: %d x %d\n", Nx, Ny);
    printf("L1 Error (Rho): %.6e\n", error_L1);
    printf("L2 Error (Rho): %.6e\n", error_L2);
    printf("======================================================\n");
}

int main(void) {
 //_control87(0, _MCW_EM); 
   // _control87(~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW), _MCW_EM);
    init_quadrature(); precompute_matrices(); init_condition();
    double t = 0.0; int nit = 0;
    printf("%-10s  %-14s\n", "Step", "Time");
    printf("--------------------------\n");

    while (t < T_END) {
        double dt = compute_dt();
        if (t + dt > T_END) dt = T_END - t;
         
        rk3_step(dt, nit);
        t += dt; nit++;
        printf("%-10d  %-14.6e %-14.6e\n", nit, t,dt);
    }
    output_results(T_END);

    compute_error(T_END); // 新增调用

    return 0;
}