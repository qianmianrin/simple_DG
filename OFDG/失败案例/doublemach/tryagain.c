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

#define N_MODE 6   
#define N_QUAD 9   

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
 * 物理安全开方函数 (保留物理底线：过大负压直接 NaN 崩溃)
 * ============================================================ */
static inline double safe_sqrt(double x) {
    if (x < -1e-10) {
        return sqrt(x); 
    }
    return sqrt(fabs(x));
}

static const double nodes1D[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

double r_quad[N_QUAD], s_quad[N_QUAD], w_quad[N_QUAD]; 

double M_diag_inv[N_MODE];          
double phi_vol[N_MODE][N_QUAD];      
double dphi_dr_vol[N_MODE][N_QUAD];  
double dphi_ds_vol[N_MODE][N_QUAD];  
double phi_face_T[N_MODE][3];        
double phi_face_B[N_MODE][3];
double phi_face_L[N_MODE][3];
double phi_face_R[N_MODE][3];

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

static const int mk_map[N_MODE] = {0, 1, 0, 2, 1, 0};
static const int nk_map[N_MODE] = {0, 0, 1, 0, 1, 2};

void precompute_matrices(void) {
    for (int k = 0; k < N_MODE; k++) {
        int mk = mk_map[k]; 
        int nk = nk_map[k]; 
        
        double m_val = (2.0 / (2.0 * mk + 1.0)) * (2.0 / (2.0 * nk + 1.0));
        M_diag_inv[k] = 1.0 / m_val;

        for (int q = 0; q < N_QUAD; q++) {
            double r = r_quad[q], s = s_quad[q];
            phi_vol[k][q]     = legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_dr_vol[k][q] = d_legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_ds_vol[k][q] = legendre_1d(mk, r) * d_legendre_1d(nk, s);
        }

        for (int p = 0; p < 3; p++) {
            double np = nodes1D[p];
            phi_face_T[k][p] = legendre_1d(mk, np) * legendre_1d(nk, 1.0);
            phi_face_B[k][p] = legendre_1d(mk, np) * legendre_1d(nk, -1.0);
            phi_face_L[k][p] = legendre_1d(mk, -1.0) * legendre_1d(nk, np);
            phi_face_R[k][p] = legendre_1d(mk, 1.0) * legendre_1d(nk, np);
        }
    }
}

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
    double c_jet = safe_sqrt(gamma * 0.4127 / 5.0);
    double max_v = 0.0;
    
    for (int i = 0; i < Nx * Ny; i++) {
        double rho = Mesh[i].U[0][0];
        double rhou = Mesh[i].U[1][0];
        double rhov = Mesh[i].U[2][0];
        double E = Mesh[i].U[3][0];
        
        double u = rhou / rho;
        double v = rhov / rho;
        double p = calc_pressure(rho, rhou, rhov, E);
        double c = safe_sqrt(gamma * p / rho);
        
        max_v = fmax(fmax(max_v, fabs(u) + c),0);
        max_v = fmax(max_v, fabs(v) + c);
    }
    return (max_v < 1e-9) ? 1.0 : max_v; 
}

static inline void llf_flux_vector(double UL[NUM_VARS], double UR[NUM_VARS], double nx, double ny, double flux_res[NUM_VARS]) {
    double FL[NUM_VARS], GL[NUM_VARS];
    double FR[NUM_VARS], GR[NUM_VARS];
    
    euler_flux(UL, FL, GL);
    euler_flux(UR, FR, GR);
    
    double rhoL = UL[0], uL = UL[1]/rhoL, vL = UL[2]/rhoL;
    double pL = calc_pressure(rhoL, UL[1], UL[2], UL[3]);
    double cL = safe_sqrt(gamma * pL / rhoL);
    double unL = uL * nx + vL * ny;
    
    double rhoR = UR[0], uR = UR[1]/rhoR, vR = UR[2]/rhoR;
    double pR = calc_pressure(rhoR, UR[1], UR[2], UR[3]);
    double cR = safe_sqrt(gamma * pR / rhoR);
    double unR = uR * nx + vR * ny;
    
    double alpha = fmax(fabs(unL) + cL, fabs(unR) + cR);
    
    for(int v = 0; v < NUM_VARS; v++) {
        double flux_n_L = FL[v] * nx + GL[v] * ny;
        double flux_n_R = FR[v] * nx + GR[v] * ny;
        flux_res[v] = 0.5 * (flux_n_L + flux_n_R - alpha * (UR[v] - UL[v]));
    }
}

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
 * 精确纠正后的特征分解矩阵 (严格遵守 L * R = I 原则)
 * ============================================================ */
static void build_eigen_x(double rho, double u, double v, double p, double L[NUM_VARS][NUM_VARS]) {
    double c = safe_sqrt(gamma * p / rho);
    if (c < 1e-14) c = 1e-14;

    double B1 = gamma - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    double inv_c2 = 1.0 / (c * c);
    
    // 行 1: 对应反向声波 (必须按 1/(2c^2) 独立缩放)
    double f_ac = 0.5 * inv_c2;
    L[0][0] = f_ac * (B2 + u * c);
    L[0][1] = f_ac * -(B1 * u + c);
    L[0][2] = f_ac * -(B1 * v);
    L[0][3] = f_ac * B1;
    
    // 行 2: 对应熵波 (必须按 1/c^2 独立缩放)
    double f_en = inv_c2;
    L[1][0] = f_en * (c * c - B2);
    L[1][1] = f_en * (B1 * u);
    L[1][2] = f_en * (B1 * v);
    L[1][3] = f_en * -B1;
    
    // 行 3: 对应剪切波 (无需除以 c^2，物理量级仅需要速度倒数 1/c 等效，或无量纲化)
    L[2][0] = -v;
    L[2][1] = 0.0;
    L[2][2] = 1.0;
    L[2][3] = 0.0;
    
    // 行 4: 对应正向声波 (必须按 1/(2c^2) 独立缩放)
    L[3][0] = f_ac * (B2 - u * c);
    L[3][1] = f_ac * -(B1 * u - c);
    L[3][2] = f_ac * -(B1 * v);
    L[3][3] = f_ac * B1;
}

static void build_eigen_y(double rho, double u, double v, double p, double L[NUM_VARS][NUM_VARS]) {
    double c = safe_sqrt(gamma * p / rho);
    if (c < 1e-14) c = 1e-14;

    double B1 = gamma - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    double inv_c2 = 1.0 / (c * c);
    
    double f_ac = 0.5 * inv_c2;
    L[0][0] = f_ac * (B2 + v * c);
    L[0][1] = f_ac * -(B1 * u);
    L[0][2] = f_ac * -(B1 * v + c);
    L[0][3] = f_ac * B1;
    
    double f_en = inv_c2;
    L[1][0] = f_en * (c * c - B2);
    L[1][1] = f_en * (B1 * u);
    L[1][2] = f_en * (B1 * v);
    L[1][3] = f_en * -B1;
    
    // 剪切波: u, -1, 0, 0
    L[2][0] = u;
    L[2][1] = -1.0;
    L[2][2] = 0.0;
    L[2][3] = 0.0;
    
    L[3][0] = f_ac * (B2 - v * c);
    L[3][1] = f_ac * -(B1 * u);
    L[3][2] = f_ac * -(B1 * v - c);
    L[3][3] = f_ac * B1;
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

    double d_dx = 2.0 / dx;
    double d_dy = 2.0 / dy;

    for (int v = 0; v < NUM_VARS; v++) {
        for (int d = 0; d < 6; d++) deriv_out[d][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            int mk = mk_map[k], nk = nk_map[k];
            double u_k = cell->U[v][k]; 
            
            deriv_out[0][v] += u_k * P_xi[mk] * P_eta[nk];
            deriv_out[1][v] += u_k * dP_xi[mk] * P_eta[nk] * d_dx;
            deriv_out[2][v] += u_k * P_xi[mk] * dP_eta[nk] * d_dy;
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

static double cell_vertex_derivs[Nx * Ny][4][6][NUM_VARS]; 
static double damp_local[Nx * Ny][3];                      

static inline void get_vertex_deriv(int ii, int jj, int vert, double out[6][NUM_VARS]) {
    if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) {
        memcpy(out, cell_vertex_derivs[jj * Nx + ii][vert], 6 * NUM_VARS * sizeof(double));
    } else {
        Element *ghost = NULL;
        if (jj < 0) ghost = &Ghost_bottom[ii];
        else if (jj >= Ny) ghost = &Ghost_top[ii];
        else if (ii < 0) ghost = &Ghost_left[jj];
        else if (ii >= Nx) ghost = &Ghost_right[jj];
        
        double xi = (vert == 1 || vert == 2) ? 1.0 : -1.0;
        double eta = (vert == 2 || vert == 3) ? 1.0 : -1.0;
        eval_element_derivatives(ghost, xi, eta, out);
    }
}

void get_element_damp_coefficients(int ii, int jj, double damp[3]) {
    Element *C_cell = &Mesh[jj * Nx + ii];
    Element *B_cell = (jj == 0) ? &Ghost_bottom[ii] : &Mesh[(jj - 1) * Nx + ii];
    Element *T_cell = (jj == Ny - 1) ? &Ghost_top[ii] : &Mesh[(jj + 1) * Nx + ii];
    Element *L_cell = (ii == 0) ? &Ghost_left[jj] : &Mesh[jj * Nx + (ii - 1)];
    Element *R_cell = (ii == Nx - 1) ? &Ghost_right[jj] : &Mesh[jj * Nx + (ii + 1)];

    double deriv_curr[4][6][NUM_VARS]; 
    double deriv_B[2][6][NUM_VARS], deriv_T[2][6][NUM_VARS];    
    double deriv_L[2][6][NUM_VARS], deriv_R[2][6][NUM_VARS];    

    get_vertex_deriv(ii, jj, 0, deriv_curr[0]); 
    get_vertex_deriv(ii, jj, 1, deriv_curr[1]); 
    get_vertex_deriv(ii, jj, 2, deriv_curr[2]); 
    get_vertex_deriv(ii, jj, 3, deriv_curr[3]); 

    get_vertex_deriv(ii, jj - 1, 3, deriv_B[0]); 
    get_vertex_deriv(ii, jj - 1, 2, deriv_B[1]); 
    
    get_vertex_deriv(ii, jj + 1, 0, deriv_T[0]); 
    get_vertex_deriv(ii, jj + 1, 1, deriv_T[1]); 
    
    get_vertex_deriv(ii - 1, jj, 1, deriv_L[0]); 
    get_vertex_deriv(ii - 1, jj, 2, deriv_L[1]); 
    
    get_vertex_deriv(ii + 1, jj, 0, deriv_R[0]); 
    get_vertex_deriv(ii + 1, jj, 3, deriv_R[1]); 

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

    for (int vert = 0; vert < 4; vert++) {
        for (int edge = 0; edge < 2; edge++) {
            double U_mean_curr[NUM_VARS], U_mean_neb[NUM_VARS];

            for (int v = 0; v < NUM_VARS; v++) {
                U_mean_curr[v] = C_cell->U[v][0];
                if (edge == 0) { 
                    if (vert == 0 || vert == 3) U_mean_neb[v] = L_cell->U[v][0];
                    else U_mean_neb[v] = R_cell->U[v][0];
                } else { 
                    if (vert == 0 || vert == 1) U_mean_neb[v] = B_cell->U[v][0];
                    else U_mean_neb[v] = T_cell->U[v][0];
                }
            }

        double rho_curr, u_curr, v_curr, p_curr;
double rho_neb, u_neb, v_neb, p_neb;

// 1. 获取两侧的原始变量
get_primitive_vars(U_mean_curr, &rho_curr, &u_curr, &v_curr, &p_curr);
get_primitive_vars(U_mean_neb, &rho_neb, &u_neb, &v_neb, &p_neb);

// 2. 计算两侧的总焓 H = (E + p) / rho
// 注意：U[3] 通常是总能 E
double H_curr = (U_mean_curr[3] + p_curr) / rho_curr;
double H_neb  = (U_mean_neb[3] + p_neb) / rho_neb;

// 3. 计算 Roe 加权因子 R = sqrt(rho_R / rho_L)
double R = sqrt(rho_neb / rho_curr);
double inv_R_plus_1 = 1.0 / (R + 1.0);

// 4. 计算 Roe 平均值
double rho_roe = R * rho_curr; // 或者 sqrt(rho_curr * rho_neb)
double u_roe   = (u_curr + R * u_neb) * inv_R_plus_1;
double v_roe   = (v_curr + R * v_neb) * inv_R_plus_1;
double H_roe   = (H_curr + R * H_neb) * inv_R_plus_1;

// 5. 从总焓 H 反推压力 p_roe (基于理想气体状态方程)
// 声速平方 c^2 = (gamma - 1) * (H - 0.5 * V^2)
double vel_sq = u_roe * u_roe + v_roe * v_roe;
double c_roe_sq = (gamma - 1.0) * (H_roe - 0.5 * vel_sq);
double p_roe = c_roe_sq * rho_roe / gamma;
   double rho_avg, u_avg, v_avg, p_avg;
// 赋值给平均变量
rho_avg = rho_roe;
u_avg   = u_roe;
v_avg   = v_roe;
p_avg   = p_roe;
    
            if (edge == 0) {
                build_eigen_x(rho_avg, u_avg, v_avg, p_avg, L_mats[vert][edge]);
            } else {
                build_eigen_y(rho_avg, u_avg, v_avg, p_avg, L_mats[vert][edge]);
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
                double alpha_norm = safe_sqrt(sum_vertices / 4.0);
                sum_alpha += alpha_norm;
            }
            if (sum_alpha > max_wave_jump) max_wave_jump = sum_alpha;
        }
        
        double factorial = 1.0;
        for (int i = 1; i <= l; i++) factorial *= i;
        double h = dx;
        double h_pow = 1.0;
        if (l == 1) h_pow = h;
        else if (l == 2) h_pow = h * h;
        
        double coeff = (2.0 * (2 * l + 1)) / (2 * k - 1) * h_pow / factorial;
        damp[l] = coeff * max_wave_jump; 
    }
}

void compute_all_damp_coefficients(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            eval_element_derivatives(&Mesh[idx], -1.0, -1.0, cell_vertex_derivs[idx][0]);
            eval_element_derivatives(&Mesh[idx],  1.0, -1.0, cell_vertex_derivs[idx][1]);
            eval_element_derivatives(&Mesh[idx],  1.0,  1.0, cell_vertex_derivs[idx][2]);
            eval_element_derivatives(&Mesh[idx], -1.0,  1.0, cell_vertex_derivs[idx][3]);
        }
    }
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            get_element_damp_coefficients(ii, jj, damp_local[idx]);
        }
    }
}

void apply_ofdg_damping(double RHS[Nx * Ny][NUM_VARS][N_MODE]) {
    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    double h_K = fmax(dx, dy); 
    
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            double *damp = damp_local[idx]; 
            
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
                llf_flux_vector(U_minus, U_p, 0, -1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_B[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
                llf_flux_vector(U_minus, U_p, 0, 1, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_T[k][p] * (dx / 2.0);
            }
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
                llf_flux_vector(U_minus, U_p, -1, 0, num_f);
                for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
                    Surf_Int[v][k] += weights1D[p] * num_f[v] * phi_face_L[k][p] * (dy / 2.0);
            }
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int p = 0; p < 3; p++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
                llf_flux_vector(U_minus, U_p, 1, 0, num_f);
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

static double compute_dt(void) { 
    return CFL * fmin(dx, dy) / (max_speed()); 
}
static inline double calc_pressure_safe(double rho, double rhou, double rhov, double E) {
    if (rho <= 0.0) return -1.0; 
    return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
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

static double U0[Nx * Ny][NUM_VARS][N_MODE], U1[Nx * Ny][NUM_VARS][N_MODE];
static double U2[Nx * Ny][NUM_VARS][N_MODE], RHS_buf[Nx * Ny][NUM_VARS][N_MODE];

/* ============================================================
 * 优化版 RK3: 冻结全局阻尼系数并彻底移除硬截断
 * ============================================================ */
void rk3_step(double dt) {
    int total = Nx * Ny;
    double h_K = fmax(dx, dy); 
    for (int i = 0; i < total; i++) memcpy(U0[i], Mesh[i].U, NUM_VARS * N_MODE * sizeof(double));
    
    // ===============================================================
    // 全局冻结 a0: 在时间步起点评定，覆盖所有子空间，绝对稳定
    // ===============================================================
    apply_ghost_cells();
    compute_all_damp_coefficients(); 
    
    double a0_global = 0.0;
    for (int i = 0; i < total; i++) {
        double val = (damp_local[i][0] + damp_local[i][1] + damp_local[i][2]) / h_K;
        if (val > a0_global) {
            a0_global = val;
        }
    }
    
    // 【Stage 1】
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 

    double z = a0_global * dt;
    double s1_const = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;

    for (int i = 0; i < total; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            U1[i][v][0] = U0[i][v][0] + dt * RHS_buf[i][v][0];
            for (int q = 1; q < N_MODE; q++) {
                U1[i][v][q] = (1.0 / s1_const) * (U0[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U0[i][v][q]);
            }
        }
    }
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U1[i], NUM_VARS * N_MODE * sizeof(double));
    apply_zhang_shu_limiter(); 
    // 【Stage 2】
    compute_rhs(RHS_buf);
    compute_all_damp_coefficients(); // 依据新流场更新阻尼
    apply_ofdg_damping(RHS_buf); 
    
    double s2_const = 1.0 + 0.5 * z + 0.125 * z * z + (1.0 / 48.0) * z * z * z;

    for (int i = 0; i < total; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            U2[i][v][0] = 0.75 * U0[i][v][0] + 0.25 * (U1[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++) {
                U2[i][v][q] = (0.75 / s2_const) * U0[i][v][q] + 
                              (s1_const / (4.0 * s2_const)) * (U1[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U1[i][v][q]);
            }
        }
    }
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U2[i], NUM_VARS * N_MODE * sizeof(double));
    apply_zhang_shu_limiter(); 
    // 【Stage 3】
    compute_rhs(RHS_buf);
    compute_all_damp_coefficients(); // 再次更新
    apply_ofdg_damping(RHS_buf); 
    
    for (int i = 0; i < total; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            Mesh[i].U[v][0] = (1.0 / 3.0) * U0[i][v][0] + (2.0 / 3.0) * (U2[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++) {
                Mesh[i].U[v][q] = (1.0 / (3.0 * s1_const)) * U0[i][v][q] + 
                                  (2.0 * s2_const / (3.0 * s1_const)) * (U2[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U2[i][v][q]);
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