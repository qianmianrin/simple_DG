#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define M_PI 3.14159265358979323846

/* ============================================================
 * 基本全局参数
 * ============================================================ */

#define Nx 100
#define Ny 50
#define xL 0.0
#define xR 1.0
#define yL -0.25
#define yR 0.25
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma (5.0 / 3.0)
#define NUM_VARS 4 


/* 时间推进参数 */
#define CFL   0.01
#define T_END 0.01
 
double max_damp=0.0;

/* ============================================================
 * 数据结构: U 中存储的是 Modal 系数
 * 从 Q2 (9个模态) 变更为 P2 (6个模态)
 * ============================================================ */
#define N_MODE 6   // P2 空间的基函数个数
#define N_QUAD 9   // 3x3 高斯积分点个数保持不变

typedef struct {
    double U[NUM_VARS][N_MODE]; 
    double face_top[NUM_VARS][3];
    double face_bottom[NUM_VARS][3];
    double face_left[NUM_VARS][3];
    double face_right[NUM_VARS][3];
} Element;

Element Mesh[Nx * Ny];
Element Ghost_bottom[Nx]; 
Element Ghost_top[Nx];    
Element Ghost_left[Ny];   
Element Ghost_right[Ny];  

/* ============================================================
 * 3 点 Gauss-Legendre 节点与权重
 * ============================================================ */
static const double nodes1D[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

double r_quad[N_QUAD], s_quad[N_QUAD], w_quad[N_QUAD]; 

/* ============================================================
 * 预计算矩阵 (Modal DG 特化)
 * ============================================================ */
double M_diag_inv[N_MODE];          
double phi_vol[N_MODE][N_QUAD];      
double dphi_dr_vol[N_MODE][N_QUAD];  
double dphi_ds_vol[N_MODE][N_QUAD];  
double phi_face_T[N_MODE][3];        
double phi_face_B[N_MODE][3];
double phi_face_L[N_MODE][3];
double phi_face_R[N_MODE][3];

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

/* 显式定义 P2 空间的基函数组合 (满足 m + n <= 2) */
static const int mk_map[N_MODE] = {0, 1, 0, 2, 1, 0};
static const int nk_map[N_MODE] = {0, 0, 1, 0, 1, 2};

void precompute_matrices(void) {
    for (int k = 0; k < N_MODE; k++) {
        int mk = mk_map[k]; // xi 方向多项式阶数
        int nk = nk_map[k]; // eta 方向多项式阶数
        
        // 正交质量矩阵对角元素
        double m_val = (2.0 / (2.0 * mk + 1.0)) * (2.0 / (2.0 * nk + 1.0));
        M_diag_inv[k] = 1.0 / m_val;

        // 体积分配置点
        for (int q = 0; q < N_QUAD; q++) {
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

/* ============================================================
 * 最大波速估算 (加入射流初始波速的 fmax 保底)
 * ============================================================ */
static inline double max_speed() { 
    double c_jet = sqrt(gamma * 0.4127 / 5.0);
    double max_v = 0.0;
    
    for (int i = 0; i < Nx * Ny; i++) {
        double rho = Mesh[i].U[0][0];
        double rhou = Mesh[i].U[1][0];
        double rhov = Mesh[i].U[2][0];
        double E = Mesh[i].U[3][0];
        
        double u = rhou / rho;
        double v = rhov / rho;
        double p = calc_pressure(rho, rhou, rhov, E);
        double c = sqrt(gamma * p / rho);
        
        max_v = fmax(fmax(max_v, fabs(u) + c), 800.0 + c_jet);
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
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];

            // 清零所有模态
            for(int v = 0; v < NUM_VARS; v++) {
                for(int k = 0; k < N_MODE; k++) {
                    cell->U[v][k] = 0.0;
                }
            }

            // 背景气体状态
            double rho = 0.5, u = 0.0, v = 0.0, p = 0.4127;
            
            double Q[NUM_VARS];
            Q[0] = rho; 
            Q[1] = rho * u; 
            Q[2] = rho * v;
            Q[3] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);

            // 只给单元均值（0阶模态）赋值
            for(int var = 0; var < NUM_VARS; var++) {
                cell->U[var][0] = Q[var];
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
            for (int k = 0; k < N_MODE; k++) {
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

/* ============================================================
 * 边界条件处理 (重点修改：现在Ghost Cell也拥有完整的模态系数 U)
 * ============================================================ */
/* ============================================================
 * 边界条件处理 (为 OFDG 探测器赋予完整的虚胞模态信息)
 * ============================================================ */
void apply_ghost_cells(void) {
    // 激波后状态 (Post-shock / Jet: State 1)
    const double rho1 = 5.0, u1 = 800.0, v1 = 0.0, p1 = 0.4127;
    const double E1 = p1 / (gamma - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);

    // 激波前状态 (Pre-shock / Ambient: State 0)
    // 注意：你提供的参考代码里此处 rho0 错写成了 5.0，应与初始场保持 0.5 一致
    const double rho0 = 0.5, u0 = 0.0, v0 = 0.0, p0 = 0.4127;
    const double E0 = p0 / (gamma - 1.0);

    // 1. 设置虚胞的模态系数 U 

    // 下边界、上边界 (Outflow / Zero-gradient) - 复制紧邻边界的内部单元
    for (int i = 0; i < Nx; i++) {
        memcpy(Ghost_bottom[i].U, Mesh[0 * Nx + i].U, NUM_VARS * N_MODE * sizeof(double));
        memcpy(Ghost_top[i].U, Mesh[(Ny - 1) * Nx + i].U, NUM_VARS * N_MODE * sizeof(double));
    }
    
    // 右边界 (Outflow)
    for (int j = 0; j < Ny; j++) {
        memcpy(Ghost_right[j].U, Mesh[j * Nx + (Nx - 1)].U, NUM_VARS * N_MODE * sizeof(double));
    }

    // 左边界 (Inflow) - 根据 y 坐标注入 Jet 或 Ambient 状态
    for (int j = 0; j < Ny; j++) {
        double y_center = yL + (j + 0.5) * dy;
        
        // 强制清零流入边界的所有高阶模态
        for(int v = 0; v < NUM_VARS; v++) {
            for(int k = 0; k < N_MODE; k++) {
                Ghost_left[j].U[v][k] = 0.0;
            }
        }

        // 将守恒量赋予均值（0阶模态）
        if (y_center >= -0.05 && y_center <= 0.05) {
            Ghost_left[j].U[0][0] = rho1;
            Ghost_left[j].U[1][0] = rho1 * u1;
            Ghost_left[j].U[2][0] = rho1 * v1;
            Ghost_left[j].U[3][0] = E1;
        } else {
            Ghost_left[j].U[0][0] = rho0;
            Ghost_left[j].U[1][0] = rho0 * u0;
            Ghost_left[j].U[2][0] = rho0 * v0;
            Ghost_left[j].U[3][0] = E0;
        }
    }

    // 2. 重新计算 Ghost 单元在各个面上的高斯点物理值（供数值通量使用）
    for (int i = 0; i < Nx; i++) {
        compute_face_values(&Ghost_bottom[i]);
        compute_face_values(&Ghost_top[i]);
    }
    for (int j = 0; j < Ny; j++) {
        compute_face_values(&Ghost_left[j]);
        compute_face_values(&Ghost_right[j]);
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
 * 激波探测器辅助工具与矩阵计算
 * ============================================================ */
static inline double calc_pressure_safe(double rho, double rhou, double rhov, double E) {
    if (rho <= 0.0) return -1.0; 
    return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
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

static void eval_legendre_basis_1d(double x, double P[3], double dP[3], double ddP[3]) {
    P[0] = 1.0; dP[0] = 0.0; ddP[0] = 0.0;
    P[1] = x;   dP[1] = 1.0; ddP[1] = 0.0;
    P[2] = 0.5 * (3.0 * x * x - 1.0); dP[2] = 3.0 * x; ddP[2] = 3.0;
}

static void eval_element_derivatives(Element *cell, double xi, double eta, double deriv_out[6][NUM_VARS]) {
    double P_xi[3], dP_xi[3], ddP_xi[3];
    double P_eta[3], dP_eta[3], ddP_eta[3];

    eval_legendre_basis_1d(xi, P_xi, dP_xi, ddP_xi);
    eval_legendre_basis_1d(eta, P_eta, dP_eta, ddP_eta);

    for (int v = 0; v < NUM_VARS; v++) {
        for (int d = 0; d < 6; d++) deriv_out[d][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            int mk = mk_map[k], nk = nk_map[k];
            double u_k = cell->U[v][k]; 
            deriv_out[0][v] += u_k * P_xi[mk] * P_eta[nk];
            deriv_out[1][v] += u_k * dP_xi[mk] * P_eta[nk];
            deriv_out[2][v] += u_k * P_xi[mk] * dP_eta[nk];
            deriv_out[3][v] += u_k * ddP_xi[mk] * P_eta[nk];
            deriv_out[4][v] += u_k * dP_xi[mk] * dP_eta[nk];
            deriv_out[5][v] += u_k * P_xi[mk] * ddP_eta[nk];
        }
    }
}

static void get_primitive_vars(double U[NUM_VARS], double *rho, double *u, double *v, double *p) {
    *rho = U[0]; *u = U[1] / U[0]; *v = U[2] / U[0];
    *p = calc_pressure(U[0], U[1], U[2], U[3]);
}

/* ============================================================
 * 核心激波探测器：计算指定单元 (ii, jj) 的 3 个 damp 系数
 * damp[0] 对应 0 阶，damp[1] 对应 1 阶，damp[2] 对应 2 阶
 * ============================================================ */
void get_element_damp_coefficients(int ii, int jj, double damp[3]) {
    Element *C_cell = &Mesh[jj * Nx + ii];
    Element *B_cell = (jj == 0) ? &Ghost_bottom[ii] : &Mesh[(jj - 1) * Nx + ii];
    Element *T_cell = (jj == Ny - 1) ? &Ghost_top[ii] : &Mesh[(jj + 1) * Nx + ii];
    Element *L_cell = (ii == 0) ? &Ghost_left[jj] : &Mesh[jj * Nx + (ii - 1)];
    Element *R_cell = (ii == Nx - 1) ? &Ghost_right[jj] : &Mesh[jj * Nx + (ii + 1)];

    double deriv_curr[4][6][NUM_VARS]; 
    double deriv_B[2][6][NUM_VARS], deriv_T[2][6][NUM_VARS];    
    double deriv_L[2][6][NUM_VARS], deriv_R[2][6][NUM_VARS];    

    eval_element_derivatives(C_cell, -1.0, -1.0, deriv_curr[0]); // v_0
    eval_element_derivatives(C_cell,  1.0, -1.0, deriv_curr[1]); // v_1
    eval_element_derivatives(C_cell,  1.0,  1.0, deriv_curr[2]); // v_2
    eval_element_derivatives(C_cell, -1.0,  1.0, deriv_curr[3]); // v_3

    eval_element_derivatives(B_cell, -1.0,  1.0, deriv_B[0]); 
    eval_element_derivatives(B_cell,  1.0,  1.0, deriv_B[1]); 
    eval_element_derivatives(T_cell, -1.0, -1.0, deriv_T[0]); 
    eval_element_derivatives(T_cell,  1.0, -1.0, deriv_T[1]); 
    eval_element_derivatives(L_cell,  1.0, -1.0, deriv_L[0]); 
    eval_element_derivatives(L_cell,  1.0,  1.0, deriv_L[1]); 
    eval_element_derivatives(R_cell, -1.0, -1.0, deriv_R[0]); 
    eval_element_derivatives(R_cell, -1.0,  1.0, deriv_R[1]); 

    double jump[4][2][6][NUM_VARS];
    for(int d=0; d<6; d++) for(int v=0; v<NUM_VARS; v++) {
        jump[0][0][d][v] = deriv_curr[0][d][v] - deriv_L[0][d][v];
        jump[0][1][d][v] = deriv_curr[0][d][v] - deriv_B[0][d][v];
        jump[1][0][d][v] = deriv_curr[1][d][v] - deriv_R[0][d][v];
        jump[1][1][d][v] = deriv_curr[1][d][v] - deriv_B[1][d][v];
        jump[2][0][d][v] = deriv_curr[2][d][v] - deriv_R[1][d][v];
        jump[2][1][d][v] = deriv_curr[2][d][v] - deriv_T[1][d][v];
        jump[3][0][d][v] = deriv_curr[3][d][v] - deriv_L[1][d][v];
        jump[3][1][d][v] = deriv_curr[3][d][v] - deriv_T[0][d][v];
    }

    double L_mats[4][2][NUM_VARS][NUM_VARS];
    double R_dummy[NUM_VARS][NUM_VARS]; 

    for (int vert = 0; vert < 4; vert++) {
        for (int edge = 0; edge < 2; edge++) {
            double U_avg[NUM_VARS];
            for (int v = 0; v < NUM_VARS; v++) {
                double u_curr = deriv_curr[vert][0][v];
                double u_neb = 0;
                if (vert == 0) u_neb = (edge == 0) ? deriv_L[0][0][v] : deriv_B[0][0][v];
                if (vert == 1) u_neb = (edge == 0) ? deriv_R[0][0][v] : deriv_B[1][0][v];
                if (vert == 2) u_neb = (edge == 0) ? deriv_R[1][0][v] : deriv_T[1][0][v];
                if (vert == 3) u_neb = (edge == 0) ? deriv_L[1][0][v] : deriv_T[0][0][v];
                U_avg[v] = 0.5 * (u_curr + u_neb);
            }
            double rho, u, v, p;
            get_primitive_vars(U_avg, &rho, &u, &v, &p);
            if (edge == 0) build_eigen_x(rho, u, v, p, R_dummy, L_mats[vert][edge]);
            else           build_eigen_y(rho, u, v, p, R_dummy, L_mats[vert][edge]);
        }
    }

    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    int k = 2; 

    for (int l = 0; l <= 2; l++) {
        double max_wave_jump = 0.0;

        for (int s = 0; s < NUM_VARS; s++) {
            double sum_alpha = 0.0; 
            for (int d = 0; d < 6; d++) {
                if (order_map[d] != l) continue;
                double sum_vertices = 0.0; 
                for (int vert = 0; vert < 4; vert++) {
                    double vertex_jump_sq = 0.0;
                    for (int edge = 0; edge < 2; edge++) {
                        double char_jump = 0.0;
                        for (int v_idx = 0; v_idx < NUM_VARS; v_idx++) {
                            char_jump += L_mats[vert][edge][s][v_idx] * jump[vert][edge][d][v_idx];
                        }
                        vertex_jump_sq += char_jump * char_jump;
                    }
                    sum_vertices += vertex_jump_sq;
                }
                double alpha_norm = sqrt(sum_vertices / 4.0);
                sum_alpha += alpha_norm;
            }
            if (sum_alpha > max_wave_jump) max_wave_jump = sum_alpha;
        }
        
        double factorial = 1.0;
        for (int i = 1; i <= l; i++) factorial *= i;
        double coeff = (2.0 * (2 * l + 1)) / (2 * k - 1) * pow(dx, l) / factorial;
        damp[l] = coeff * max_wave_jump; 
    }
}
/* ============================================================
 * 嵌合 OFDG 阻尼项：根据探测器得出的 damp 系数，衰减高阶 RHS
 * ============================================================ */
void apply_ofdg_damping(double RHS[Nx * Ny][NUM_VARS][N_MODE]) {
    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            
            // 1. 调用你原本代码中的探测器，获取该单元 0, 1, 2 阶的阻尼系数
            get_element_damp_coefficients(ii, jj, damp);
            
            // 2. 物理守恒性：强制 0 阶（单元平均值）的阻尼为 0
            damp[0] = 0.0; 
double temp=damp[1]+damp[2]+damp[0];
max_damp=fmax(temp,max_damp);
            // 3. 对 1 阶和 2 阶的模态 RHS 进行阻尼修正
            for (int v = 0; v < NUM_VARS; v++) {
                for (int k = 1; k < N_MODE; k++) { // k 从 1 开始，保留 k=0 的守恒量
                    int l = order_map[k];
                    double h=dx;
                    RHS[idx][v][k] -=(damp[l]/h) * Mesh[idx].U[v][k];
                }
            }
        }
    }
}
/* ============================================================
 * 正保真度限制器 (Positivity Limiter)
 * ============================================================ */
void apply_positivity_limiter(void) {
    const double eps = 1e-13;

    for (int idx = 0; idx < Nx * Ny; idx++) {
        Element *cell = &Mesh[idx];

        double Ubar[NUM_VARS];
        for (int v = 0; v < NUM_VARS; v++) Ubar[v] = cell->U[v][0];

        if (Ubar[0] < eps) {
            for (int k = 1; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++) cell->U[v][k] = 0.0;
            continue;
        }

        double U_test[21][NUM_VARS];
        int pt = 0;

        for (int q = 0; q < N_QUAD; q++, pt++) {
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_vol[k][q];
        }
        for (int q = 0; q < 3; q++) {
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_T[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_B[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_L[k][q];
            pt++;
            for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[pt][v] += cell->U[v][k] * phi_face_R[k][q];
            pt++;
        }

        double rho_min = U_test[0][0];
        for (int i = 1; i < 21; i++)
            if (U_test[i][0] < rho_min) rho_min = U_test[i][0];

        double theta1 = 1.0;
        if (rho_min < eps) {
            theta1 = (Ubar[0] - eps) / (Ubar[0] - rho_min);
            if (theta1 < 0.0) theta1 = 0.0;
            if (theta1 > 1.0) theta1 = 1.0;

            for (int v = 0; v < NUM_VARS; v++)
                for (int k = 1; k < N_MODE; k++)
                    cell->U[v][k] *= theta1;

            for (int i = 0; i < 21; i++)
                for (int v = 0; v < NUM_VARS; v++)
                    U_test[i][v] = Ubar[v] + theta1 * (U_test[i][v] - Ubar[v]);
        }

        double theta2 = 1.0;
        for (int i = 0; i < 21; i++) {
            double p_i = calc_pressure_safe(U_test[i][0], U_test[i][1],
                                            U_test[i][2], U_test[i][3]);
            if (p_i < eps) {
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
                double t_star = t_L; 
                if (t_star < theta2) theta2 = t_star;
            }
        }

        if (theta2 < 1.0)
            for (int v = 0; v < NUM_VARS; v++)
                for (int k = 1; k < N_MODE; k++)
                    cell->U[v][k] *= theta2;
    }
}

/* ============================================================
 * 弱形式 DG 右端项 (Weak Form)
 * ============================================================ */
void compute_rhs(double RHS[Nx * Ny][NUM_VARS][N_MODE]) {
    boundary_value();
    apply_ghost_cells();

    double J = (dx * dy) / 4.0; 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double u_plus[NUM_VARS][3];

            double Vol_Int[NUM_VARS][N_MODE] = {{0}};
            double Surf_Int[NUM_VARS][N_MODE] = {{0}};

            for (int q = 0; q < N_QUAD; q++) {
                double U_phys[NUM_VARS] = {0};
                for (int k = 0; k < N_MODE; k++)
                    for(int v=0; v<NUM_VARS; v++) U_phys[v] += cell->U[v][k] * phi_vol[k][q];
                
                double F_val[NUM_VARS], G_val[NUM_VARS];
                euler_flux(U_phys, F_val, G_val);

                for (int k = 0; k < N_MODE; k++) {
                    for(int v=0; v<NUM_VARS; v++) {
                        Vol_Int[v][k] += w_quad[q] * ( F_val[v] * dphi_dr_vol[k][q] * (dy / 2.0) + 
                                                       G_val[v] * dphi_ds_vol[k][q] * (dx / 2.0) );
                    }
                }
            }

            double U_minus[NUM_VARS], U_p[NUM_VARS], num_f[NUM_VARS];

            get_neighbor_face(ii, jj, 0, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_bottom[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 0, -1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_B[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 0, 1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_T[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, -1, 0, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_L[k][p] * (dy / 2.0);
            }
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
                lf_flux_vector(U_minus, U_p, 1, 0, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_R[k][p] * (dy / 2.0);
            }

            for(int v=0; v<NUM_VARS; v++) {
                for (int k = 0; k < N_MODE; k++) {
                    RHS[idx][v][k] = (Vol_Int[v][k] - Surf_Int[v][k]) * M_diag_inv[k] / J;
                }
            }
        }
    }
}

/* ============================================================
 * 时间推进 
 * ============================================================ */
static double compute_dt(void) { 
    double deltat= CFL * fmin(dx, dy) / (max_speed()+max_damp); 
max_damp=0;
return deltat;
}

static double U0[Nx * Ny][NUM_VARS][N_MODE], U1[Nx * Ny][NUM_VARS][N_MODE];
static double U2[Nx * Ny][NUM_VARS][N_MODE], RHS_buf[Nx * Ny][NUM_VARS][N_MODE];

/* ============================================================
 * 完整嵌合了 OFDG 指数型 TVD-RK3 的时间推进
 * ============================================================ */
/* ============================================================
 * 完整嵌合了 OFDG 修正指数型 TVD-RK3 的时间推进 (严格对应论文 2.24 格式)
 * ============================================================ */
void rk3_step(double dt) {
    int total = Nx * Ny;
    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    
    // 备份 U^n
    for (int i = 0; i < total; i++) {
        memcpy(U0[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));
    }

    // -----------------------------------------------------------------
    // Stage 1
    compute_rhs(RHS_buf); 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            get_element_damp_coefficients(ii, jj, damp);
            damp[0] = 0.0; // 0阶模态（平均值）绝不衰减，保证守恒

            for (int v = 0; v < NUM_VARS; v++) {
                for (int k = 0; k < N_MODE; k++) {
                    double z = damp[order_map[k]] * dt;
                    double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
                    
                    // 对应论文 (2.24a)
                    U1[idx][v][k] = (U0[idx][v][k] + dt * RHS_buf[idx][v][k]) / s1;
                }
            }
        }
    }
    
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U1[i], NUM_VARS * N_MODE * sizeof(double));
    apply_positivity_limiter(); 
    // 必须将 limiter 截断后的安全值拷贝回 U1 供下一阶段使用
    for (int i = 0; i < total; i++) memcpy(U1[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double)); 

    // -----------------------------------------------------------------
    // Stage 2
    compute_rhs(RHS_buf); 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            get_element_damp_coefficients(ii, jj, damp);
            damp[0] = 0.0;

            for (int v = 0; v < NUM_VARS; v++) {
                for (int k = 0; k < N_MODE; k++) {
                    double z = damp[order_map[k]] * dt;
                    double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
                    double s2 = 1.0 + 0.5 * z + (1.0 / 8.0) * z * z + (1.0 / 48.0) * z * z * z;
                    
                    // 对应论文 (2.24b)
                    double term1 = (0.75 / s2) * U0[idx][v][k];
                    double term2 = (s1 / (4.0 * s2)) * (U1[idx][v][k] + dt * RHS_buf[idx][v][k]);
                    
                    U2[idx][v][k] = term1 + term2;
                }
            }
        }
    }
    
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U2[i], NUM_VARS * N_MODE * sizeof(double));
   apply_positivity_limiter(); 
    for (int i = 0; i < total; i++) memcpy(U2[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));

    // -----------------------------------------------------------------
    // Stage 3
    compute_rhs(RHS_buf); 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            get_element_damp_coefficients(ii, jj, damp);
            damp[0] = 0.0;

            for (int v = 0; v < NUM_VARS; v++) {
                for (int k = 0; k < N_MODE; k++) {
                    double z = damp[order_map[k]] * dt;
                    double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
                    double s2 = 1.0 + 0.5 * z + (1.0 / 8.0) * z * z + (1.0 / 48.0) * z * z * z;
                    
                    // 对应论文 (2.24c)
                    double term1 = (1.0 / (3.0 * s1)) * U0[idx][v][k];
                    double term2 = (2.0 * s2 / (3.0 * s1)) * (U2[idx][v][k] + dt * RHS_buf[idx][v][k]);
                    
                    Mesh[idx].U[v][k] = term1 + term2;
                }
            }
        }
    }
    
   apply_positivity_limiter(); 
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

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double r = nodes1D[j], s = nodes1D[i];
                    double x_phys = xc + (dx / 2.0) * r, y_phys = yc + (dy / 2.0) * s;
                    
                    double U_phys[NUM_VARS] = {0};
                    for (int k = 0; k < N_MODE; k++) {
                        double phi_val = legendre_1d(mk_map[k], r) * legendre_1d(nk_map[k], s);
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
}



int main(void) {

    init_quadrature(); precompute_matrices(); init_condition();
    double t = 0.0; int nit = 0;
    printf("%-10s  %-14s\n", "Step", "Time");
    printf("--------------------------\n");

    while (t < T_END) {
        double dt = compute_dt();
        if (t + dt > T_END) dt = T_END - t;
        rk3_step(dt); t += dt; nit++;
           printf("%-10d  %-14.6e %-14.6e %-14.6e\n", nit, t,dt,max_damp);
    }
    output_results(T_END); 
    return 0;
}