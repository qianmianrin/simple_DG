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
#define Nx 200
#define Ny 50
#define xL 0.0
#define xR 4.0
#define yL 0.0
#define yR 1.0
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma 1.4
#define NUM_VARS 4 

/* 时间推进参数 */
#define CFL   0.01
#define T_END 0.2
//双马赫反射需要全局时间
double current_time = 0.0;

double max_damp=1.0;

/* ============================================================
 * 数据结构: U 中存储的是 Modal 系数
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

static inline double max_speed() { 
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
        
        max_v = fmax(max_v, fabs(u) + c);
        max_v = fmax(max_v, fabs(v) + c);
    }
    return (max_v < 1e-9) ? 1.0 : max_v; 
}

static inline void hllc_flux_vector(double UL[NUM_VARS], double UR[NUM_VARS], double nx, double ny, double flux_res[NUM_VARS]) {
    // 1. 提取左右状态的基础物理量
    double rhoL = UL[0], uL = UL[1] / rhoL, vL = UL[2] / rhoL, EL = UL[3];
    double pL = calc_pressure(rhoL, UL[1], UL[2], EL);
    double cL = sqrt(gamma * pL / rhoL);
    double unL = uL * nx + vL * ny;
    double HL = (EL + pL) / rhoL;

    double rhoR = UR[0], uR = UR[1] / rhoR, vR = UR[2] / rhoR, ER = UR[3];
    double pR = calc_pressure(rhoR, UR[1], UR[2], ER);
    double cR = sqrt(gamma * pR / rhoR);
    double unR = uR * nx + vR * ny;
    double HR = (ER + pR) / rhoR;

    // 2. 计算基于法向的物理通量 F_n
    double F_L[NUM_VARS] = {
        rhoL * unL,
        rhoL * uL * unL + pL * nx,
        rhoL * vL * unL + pL * ny,
        unL * (EL + pL)
    };

    double F_R[NUM_VARS] = {
        rhoR * unR,
        rhoR * uR * unR + pR * nx,
        rhoR * vR * unR + pR * ny,
        unR * (ER + pR)
    };

    // 3. Roe 平均估计波速
    double R_roe = sqrt(rhoR / rhoL);
    double un_roe = (unL + R_roe * unR) / (1.0 + R_roe);
    double u_roe = (uL + R_roe * uR) / (1.0 + R_roe);
    double v_roe = (vL + R_roe * vR) / (1.0 + R_roe);
    double H_roe = (HL + R_roe * HR) / (1.0 + R_roe);
    double c2_roe = (gamma - 1.0) * (H_roe - 0.5 * (u_roe * u_roe + v_roe * v_roe));
    double c_roe = sqrt(fmax(c2_roe, 1e-12));

    // 4. 计算三种特征波速 (SL, SR, S*)
    double SL = fmin(unL - cL, un_roe - c_roe);
    double SR = fmax(unR + cR, un_roe + c_roe);

    // 中间接触间断波速 S*
    double num = pR - pL + rhoL * unL * (SL - unL) - rhoR * unR * (SR - unR);
    double den = rhoL * (SL - unL) - rhoR * (SR - unR);
    double S_star = num / den;

    // 5. 根据 HLLC 判断逻辑计算最终的交界面通量
    if (SL >= 0.0) {
        for (int v = 0; v < NUM_VARS; v++) flux_res[v] = F_L[v];
    } 
    else if (SR <= 0.0) {
        for (int v = 0; v < NUM_VARS; v++) flux_res[v] = F_R[v];
    } 
    else if (SL <= 0.0 && S_star >= 0.0) {
        // 左侧星区状态 U*L
        double inv_star_L = 1.0 / (SL - S_star);
        double U_star_L[NUM_VARS] = {
            rhoL * (SL - unL) * inv_star_L,
            rhoL * (SL - unL) * inv_star_L * (uL + (S_star - unL) * nx),
            rhoL * (SL - unL) * inv_star_L * (vL + (S_star - unL) * ny),
            rhoL * (SL - unL) * inv_star_L * (EL / rhoL + (S_star - unL) * (S_star + pL / (rhoL * (SL - unL))))
        };
        for (int v = 0; v < NUM_VARS; v++) 
            flux_res[v] = F_L[v] + SL * (U_star_L[v] - UL[v]);
    } 
    else { // S_star <= 0.0 && SR >= 0.0
        // 右侧星区状态 U*R
        double inv_star_R = 1.0 / (SR - S_star);
        double U_star_R[NUM_VARS] = {
            rhoR * (SR - unR) * inv_star_R,
            rhoR * (SR - unR) * inv_star_R * (uR + (S_star - unR) * nx),
            rhoR * (SR - unR) * inv_star_R * (vR + (S_star - unR) * ny),
            rhoR * (SR - unR) * inv_star_R * (ER / rhoR + (S_star - unR) * (S_star + pR / (rhoR * (SR - unR))))
        };
        for (int v = 0; v < NUM_VARS; v++) 
            flux_res[v] = F_R[v] + SR * (U_star_R[v] - UR[v]);
    }
}
/* ============================================================
 * 初始条件: L2 投影到 Modal 空间
 * ============================================================ */


void init_condition(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];
            double x_center = xL + (ii + 0.5) * dx;
            double y_center = yL + (jj + 0.5) * dy;

            for(int v=0; v<NUM_VARS; v++)
                for(int k=0; k<N_MODE; k++) cell->U[v][k] = 0.0;

           
                double x = x_center ;
                double y = y_center ;

                   
    double theta=M_PI/3.0;
            double shock_X=1.0/6.0+y/(tan(theta));
                double rho, u, v, p;
        
          if (x <shock_X)      { rho = 8.0;    u = 4.125*sqrt(3.0);   v = -4.125;   p = 116.5;   }
           
            else                            { rho = 1.4; u = 0.0;   v = 0.0; p = 1.0;   }
            
                double Q[NUM_VARS];
                Q[0] = rho; Q[1] = rho * u; Q[2] = rho * v;
                Q[3] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);

           
                    for(int var=0; var<NUM_VARS; var++)
                    {cell->U[var][0] = Q[var] ;}// 
            
            
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

// 边界条件 
// 
void apply_ghost_cells(void) {
    // 激波后状态 (Post-shock: State 1)
    const double rho1 = 8.0;
    const double u1   = 4.125 * sqrt(3.0);
    const double v1   = -4.125;
    const double p1   = 116.5;
    const double E1   = p1 / (gamma - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);
    const double Q1[4] = {rho1, rho1 * u1, rho1 * v1, E1};

    // 激波前状态 (Pre-shock: State 0)
    const double rho0 = 1.4;
    const double u0   = 0.0;
    const double v0   = 0.0;
    const double p0   = 1.0;
    const double E0   = p0 / (gamma - 1.0);
    const double Q0[4] = {rho0, rho0 * u0, rho0 * v0, E0};

    // 初始化：清空所有虚单元的模态系数
    memset(Ghost_bottom, 0, sizeof(Element) * Nx);
    memset(Ghost_top, 0, sizeof(Element) * Nx);
    memset(Ghost_left, 0, sizeof(Element) * Ny);
    memset(Ghost_right, 0, sizeof(Element) * Ny);

    // 1. 左边界 (Left, x = 0.0): 始终为激波后状态 (常数分布，仅 0 阶模态有值)
    for (int j = 0; j < Ny; j++) {
        for (int v = 0; v < NUM_VARS; v++) Ghost_left[j].U[v][0] = Q1[v];
    }

    // 2. 右边界 (Right, x = 1.0): 零梯度流出 (直接复制最右侧单元的模态系数)
    for (int j = 0; j < Ny; j++) {
        memcpy(Ghost_right[j].U, Mesh[j * Nx + (Nx - 1)].U, sizeof(double) * NUM_VARS * N_MODE);
    }

    // 3. 下边界 (Bottom, y = 0)
    for (int i = 0; i < Nx; i++) {
        double xc = xL + (i + 0.5) * dx;
        if (xc < 1.0 / 6.0) {
            // 激波流入区：常数状态 Q1
            for (int v = 0; v < NUM_VARS; v++) Ghost_bottom[i].U[v][0] = Q1[v];
        } else {
            // 反射边界：根据内部第一个单元进行镜像映射
            // 密度(0)、X动量(1)、能量(3)是偶对称，Y动量(2)是奇对称
            for (int k = 0; k < N_MODE; k++) {
                int nk = nk_map[k]; // 获取 y 方向的阶数
                double sign = (nk % 2 == 0) ? 1.0 : -1.0; // 考虑基函数在镜像位置的符号
                
                Ghost_bottom[i].U[0][k] =  sign * Mesh[0 * Nx + i].U[0][k];
                Ghost_bottom[i].U[1][k] =  sign * Mesh[0 * Nx + i].U[1][k];
                Ghost_bottom[i].U[2][k] = -sign * Mesh[0 * Nx + i].U[2][k]; // Y 动量取反
                Ghost_bottom[i].U[3][k] =  sign * Mesh[0 * Nx + i].U[3][k];
            }
        }
    }

    // 4. 上边界 (Top, y = 1.0): 随时间移动的精确激波解
    double sin60 = sqrt(3.0) / 2.0;
        double xs_top = (1.0 / 6.0) + (10.0 * current_time + 0.5) / sin60;
    // 注意：这里的 xs_top 逻辑需根据双马赫反射的标准定义调整，此处设为常数状态切换
    for (int i = 0; i < Nx; i++) {
        double xc = xL + (i + 0.5) * dx;
        if (xc < xs_top) {
            for (int v = 0; v < NUM_VARS; v++) Ghost_top[i].U[v][0] = Q1[v];
        } else {
            for (int v = 0; v < NUM_VARS; v++) Ghost_top[i].U[v][0] = Q0[v];
        }
    }

    // 最后统一更新所有虚单元的面值（用于计算通量）
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
    double c = sqrt(gamma * p / rho);
    double H = (p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v) + p) / rho;
    
    double B1 = gamma - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    
    // X方向法向: nx = 1, ny = 0 -> u_hat = u, v_hat = v
    double factor_L = B1 / c; 
    
    // 严格匹配论文公式 (2.14) 的左特征矩阵 L (即 R^{-1})
    L[0][0] = factor_L * 0.5 * (B2 + u * c);
    L[0][1] = factor_L * -0.5 * (B1 * u + c);
    L[0][2] = factor_L * -0.5 * (B1 * v);
    L[0][3] = factor_L * 0.5 * B1;
    
    L[1][0] = factor_L * (c * c - B2);
    L[1][1] = factor_L * B1 * u;
    L[1][2] = factor_L * B1 * v;
    L[1][3] = factor_L * -B1;
    
    L[2][0] = factor_L * (-v * c);
    L[2][1] = 0.0;
    L[2][2] = factor_L * c;
    L[2][3] = 0.0;
    
    L[3][0] = factor_L * 0.5 * (B2 - u * c);
    L[3][1] = factor_L * -0.5 * (B1 * u - c);
    L[3][2] = factor_L * -0.5 * (B1 * v);
    L[3][3] = factor_L * 0.5 * B1;

    // 构造对应的严格逆矩阵 R (即右特征矩阵，确保 L * R = I)
    double fR_c = 1.0 / (c * B1);
    double fR_1 = 1.0 / B1;
    
    R[0][0] = fR_c;            R[0][1] = fR_c;          R[0][2] = 0.0;       R[0][3] = fR_c;
    R[1][0] = fR_c * (u - c);  R[1][1] = fR_c * u;      R[1][2] = 0.0;       R[1][3] = fR_c * (u + c);
    R[2][0] = fR_c * v;        R[2][1] = fR_c * v;      R[2][2] = fR_1;      R[2][3] = fR_c * v;
    R[3][0] = fR_c * (H - u*c);R[3][1] = fR_c *(B2/B1); R[3][2] = fR_1 * v;  R[3][3] = fR_c * (H + u*c);
}

static void build_eigen_y(double rho, double u, double v, double p, double R[NUM_VARS][NUM_VARS], double L[NUM_VARS][NUM_VARS]) {
    double c = sqrt(gamma * p / rho);
    double H = (p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v) + p) / rho;
    
    double B1 = gamma - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    
    // Y方向法向: nx = 0, ny = 1 -> u_hat = v, v_hat = -u
    double factor_L = B1 / c;
    
    // 严格匹配论文公式 (2.14) 的左特征矩阵 L (即 R^{-1})
    L[0][0] = factor_L * 0.5 * (B2 + v * c);
    L[0][1] = factor_L * -0.5 * (B1 * u);
    L[0][2] = factor_L * -0.5 * (B1 * v + c);
    L[0][3] = factor_L * 0.5 * B1;
    
    L[1][0] = factor_L * (c * c - B2);
    L[1][1] = factor_L * B1 * u;
    L[1][2] = factor_L * B1 * v;
    L[1][3] = factor_L * -B1;
    
    L[2][0] = factor_L * (u * c);
    L[2][1] = factor_L * (-c);
    L[2][2] = 0.0;
    L[2][3] = 0.0;
    
    L[3][0] = factor_L * 0.5 * (B2 - v * c);
    L[3][1] = factor_L * -0.5 * (B1 * u);
    L[3][2] = factor_L * -0.5 * (B1 * v - c);
    L[3][3] = factor_L * 0.5 * B1;

    // 构造对应的严格逆矩阵 R (即右特征矩阵，确保 L * R = I)
    double fR_c = 1.0 / (c * B1);
    double fR_1 = 1.0 / B1;
    
    R[0][0] = fR_c;            R[0][1] = fR_c;          R[0][2] = 0.0;       R[0][3] = fR_c;
    R[1][0] = fR_c * u;        R[1][1] = fR_c * u;      R[1][2] = -fR_1;     R[1][3] = fR_c * u;
    R[2][0] = fR_c * (v - c);  R[2][1] = fR_c * v;      R[2][2] = 0.0;       R[2][3] = fR_c * (v + c);
    R[3][0] = fR_c * (H - v*c);R[3][1] = fR_c *(B2/B1); R[3][2] = -fR_1 * u; R[3][3] = fR_c * (H + v*c);
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

    // 引入从参考系到物理系的链式法则缩放项
    double d_dx = 2.0 / dx;
    double d_dy = 2.0 / dy;

    for (int v = 0; v < NUM_VARS; v++) {
        for (int d = 0; d < 6; d++) deriv_out[d][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            int mk = mk_map[k], nk = nk_map[k];
            double u_k = cell->U[v][k]; 
            
            // 0 阶导数 (函数值本身不受坐标缩放影响)
            deriv_out[0][v] += u_k * P_xi[mk] * P_eta[nk];
            // 1 阶导数
            deriv_out[1][v] += u_k * dP_xi[mk] * P_eta[nk] * d_dx;
            deriv_out[2][v] += u_k * P_xi[mk] * dP_eta[nk] * d_dy;
            // 2 阶导数
            deriv_out[3][v] += u_k * ddP_xi[mk] * P_eta[nk] * (d_dx * d_dx);
            deriv_out[4][v] += u_k * dP_xi[mk] * dP_eta[nk] * (d_dx * d_dy);
            deriv_out[5][v] += u_k * P_xi[mk] * ddP_eta[nk] * (d_dy * d_dy);
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
            
            double U_mean_curr[NUM_VARS], U_mean_neb[NUM_VARS];

            for (int v = 0; v < NUM_VARS; v++) {
                U_mean_curr[v] = C_cell->U[v][0];
                
                if (edge == 0) { // X 方向交界面
                    if (vert == 0 || vert == 3) { // 顶点的左侧边
                        U_mean_neb[v] = L_cell->U[v][0];
                    } else { // 顶点的右侧边
                        U_mean_neb[v] = R_cell->U[v][0];
                    }
                } else { // Y 方向交界面
                    if (vert == 0 || vert == 1) { // 顶点的下侧边
                        U_mean_neb[v] = B_cell->U[v][0];
                    } else { // 顶点的上侧边
                        U_mean_neb[v] = T_cell->U[v][0];
                    }
                }
            }
        
            
            double rhoL, uL, vL, pL, HL;
            double rhoR, uR, vR, pR, HR;
            
            get_primitive_vars(U_mean_curr, &rhoL, &uL, &vL, &pL);
            get_primitive_vars(U_mean_neb, &rhoR, &uR, &vR, &pR);
        
            HL = (U_mean_curr[3] + pL) / rhoL;
            HR = (U_mean_neb[3] + pR) / rhoR;
            
            double sqL = sqrt(rhoL);
            double sqR = sqrt(rhoR);
            double inv_sq = 1.0 / (sqL + sqR);
            
            double rho_roe = sqL * sqR;
            double u_roe   = (sqL * uL + sqR * uR) * inv_sq;
            double v_roe   = (sqL * vL + sqR * vR) * inv_sq;
            double H_roe   = (sqL * HL + sqR * HR) * inv_sq;
            
            double q2_roe = 0.5 * (u_roe * u_roe + v_roe * v_roe);
            double c2_roe = (gamma - 1.0) * (H_roe - q2_roe);
            double p_roe = rho_roe * c2_roe / gamma;
            
            // 构建标准界面特征矩阵
            if (edge == 0) {
                build_eigen_x(rho_roe, u_roe, v_roe, p_roe, R_dummy, L_mats[vert][edge]);
            } else {
                build_eigen_y(rho_roe, u_roe, v_roe, p_roe, R_dummy, L_mats[vert][edge]);
            }
            
  
           
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
        double h = dx;
       // 显式计算整数次幂，避免调用缓慢的 pow() 函数
        double h_pow = 1.0;
        if (l == 1) h_pow = h;
        else if (l == 2) h_pow = h * h;
        
        double coeff = (2.0 * (2 * l + 1)) / (2 * k - 1) * h_pow / factorial;
        damp[l] = coeff * max_wave_jump; 

   
    }

  

}
/* ============================================================
 * 嵌合 OFDG 阻尼项：根据探测器得出的 damp 系数，衰减高阶 RHS
 * ============================================================ */
void apply_ofdg_damping(double RHS[Nx * Ny][NUM_VARS][N_MODE]) {
    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    double h_K = fmax(dx, dy); 
    
    max_damp = 0.0; 
    
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            
            get_element_damp_coefficients(ii, jj, damp);
            
            // 【关键数学纠错】必须在这里除以 h_K，确保 a_0 与 Sigma 矩阵尺度严格一致
            double current_total_damp = (damp[0] + damp[1] + damp[2]) / h_K;
            if (current_total_damp > max_damp) {
                max_damp = current_total_damp;
            }
            
            double sigma_k[N_MODE];
            sigma_k[0] = 0.0; 
            for (int k = 1; k < N_MODE; k++) {
                int dk = order_map[k];
                double sum_sigma = 0.0;
                for (int l = 0; l <= dk; l++) {
                    sum_sigma += damp[l] / h_K;
                }
                sigma_k[k] = sum_sigma;
            }
            
            for (int v = 0; v < NUM_VARS; v++) {
                for (int k = 1; k < N_MODE; k++) { 
                    RHS[idx][v][k] -= sigma_k[k] * Mesh[idx].U[v][k];
                }
            }
        }
    }
}

static double legendre_eval(int order, double x) {
    if (order == 0) return 1.0;
    if (order == 1) return x;
    if (order == 2) return 0.5 * (3.0 * x * x - 1.0);
    return 0.0; 
}
void apply_zhang_shu_limiter(void) {

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
                hllc_flux_vector(U_minus, U_p, 0, -1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_B[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
                hllc_flux_vector(U_minus, U_p, 0, 1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_T[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
                hllc_flux_vector(U_minus, U_p, -1, 0, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_L[k][p] * (dy / 2.0);
            }
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
                hllc_flux_vector(U_minus, U_p, 1, 0, num_f);
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
      
    return CFL * fmin(dx, dy) / (max_speed()); }

static double U0[Nx * Ny][NUM_VARS][N_MODE], U1[Nx * Ny][NUM_VARS][N_MODE];
static double U2[Nx * Ny][NUM_VARS][N_MODE], RHS_buf[Nx * Ny][NUM_VARS][N_MODE];
static double a0_local[Nx * Ny];

/* ============================================================
 * 时间推进
 * ============================================================ */

/* 在全局变量区确保有这个声明（如果在上一轮已经加了就保留） */
static double a0_local[Nx * Ny];

/* ============================================================
 * 修正版 RK3: 局部 a0 + 严格同步的保正限制器 (Positivity-Preserving)
 * ============================================================ */
void rk3_step(double dt) {
    int total = Nx * Ny;
    for (int i = 0; i < total; i++) memcpy(U0[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));
    
    // 【阶段 0】：获取基础 RHS，同时提取每个单元的局部 a_0
    compute_rhs(RHS_buf);
    
    double h_K = fmax(dx, dy); 
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double damp[3];
            get_element_damp_coefficients(ii, jj, damp);
            // 将当前单元的最大阻尼系数作为该单元专属的局部积分因子 a0
            a0_local[idx] = (damp[0] + damp[1] + damp[2]) / h_K;
        }
    }
    
    // 应用 OFDG 阻尼到 RHS 
    apply_ofdg_damping(RHS_buf); 

    // -----------------------------------------------------------------
    // Stage 1
    for (int i = 0; i < total; i++) {
        double a0 = a0_local[i];
        double z  = a0 * dt;
        double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;

        for (int v = 0; v < NUM_VARS; v++) {
            U1[i][v][0] = U0[i][v][0] + dt * RHS_buf[i][v][0];
            for (int q = 1; q < N_MODE; q++) {
                U1[i][v][q] = (1.0 / s1) * (U0[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0 * U0[i][v][q]);
            }
        }
    }
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U1[i], NUM_VARS * N_MODE * sizeof(double));
    
    // 必须恢复保正限制器拦截 p < 0
    apply_zhang_shu_limiter(); 
    
    // 【关键修复】：限制器修改了 Mesh.U，必须同步回 U1，否则 Stage 2 仍会混入危险值
    for (int i = 0; i < total; i++) memcpy(U1[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));

    // -----------------------------------------------------------------
    // Stage 2
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 
    
    for (int i = 0; i < total; i++) {
        double a0 = a0_local[i];
        double z  = a0 * dt;
        double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
        double s2 = 1.0 + 0.5 * z + 0.125 * z * z + (1.0 / 48.0) * z * z * z;

        for (int v = 0; v < NUM_VARS; v++) {
            U2[i][v][0] = 0.75 * U0[i][v][0] + 0.25 * (U1[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++) {
                U2[i][v][q] = (0.75 / s2) * U0[i][v][q] + 
                              (s1 / (4.0 * s2)) * (U1[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0 * U1[i][v][q]);
            }
        }
    }
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U2[i], NUM_VARS * N_MODE * sizeof(double));
    
    apply_zhang_shu_limiter();
    // 【关键修复】：同步回 U2
    for (int i = 0; i < total; i++) memcpy(U2[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));

    // -----------------------------------------------------------------
    // Stage 3
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 
    
    for (int i = 0; i < total; i++) {
        double a0 = a0_local[i];
        double z  = a0 * dt;
        double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
        double s2 = 1.0 + 0.5 * z + 0.125 * z * z + (1.0 / 48.0) * z * z * z;

        for (int v = 0; v < NUM_VARS; v++) {
            Mesh[i].U[v][0] = (1.0 / 3.0) * U0[i][v][0] + (2.0 / 3.0) * (U2[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++) {
                Mesh[i].U[v][q] = (1.0 / (3.0 * s1)) * U0[i][v][q] + 
                                  (2.0 * s2 / (3.0 * s1)) * (U2[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0 * U2[i][v][q]);
            }
        }
    }
    apply_zhang_shu_limiter();
}
void output_results(double t) {
    FILE *fp = fopen("result2.dat", "w");
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
 //  _control87(0, _MCW_EM); 
   // _control87(~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW), _MCW_EM);
    init_quadrature(); precompute_matrices(); init_condition();
    double t = 0.0; int nit = 0;
    printf("%-10s  %-14s\n", "Step", "Time");
    printf("--------------------------\n");

    while (t < T_END) {
        double dt = compute_dt();
        if (t + dt > T_END) dt = T_END - t;

        current_time=t;

        rk3_step(dt); t += dt; nit++;
           printf("%-10d  %-14.6e %-14.6e \n", nit, t,dt);
      
    }
    output_results(T_END);
    return 0;
}