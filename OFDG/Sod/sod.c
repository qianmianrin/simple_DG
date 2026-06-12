#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define M_PI 3.14159265358979323846

/* ============================================================
 * 全局参数 (1D Sedov Blast Wave)
 * ============================================================ */
#define Nx 513
#define xL -2.0
#define xR 2.0
#define dx ((xR - xL) / Nx)
#define gamma 1.4
#define NUM_VARS 3   
#define N_MODE 3     
#define N_QUAD 3     

#define CFL   0.1
#define T_END 1.3

/* ============================================================
 * 数据结构
 * ============================================================ */
typedef struct {
    double U[NUM_VARS][N_MODE]; // 模态系数
    double face_L[NUM_VARS];    // 左侧面重构物理值
    double face_R[NUM_VARS];    // 右侧面重构物理值
} Element;

Element Mesh[Nx];
Element Ghost_left;
Element Ghost_right;

/* 3 点 Gauss-Legendre 节点与权重 */
static const double nodes1D[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

double M_diag_inv[N_MODE];
double phi_vol[N_MODE][N_QUAD];
double dphi_vol[N_MODE][N_QUAD];
double phi_face_L[N_MODE];
double phi_face_R[N_MODE];

/* 全局缓存 */
static double damp_local[Nx][3];

/* ============================================================
 * 物理安全开方函数
 * ============================================================ */
static inline double safe_sqrt(double x) {
    if (x < -1e-10) {
        // 超过微小机器误差的负数，直接调用 sqrt 产生 NaN，让程序失去物理意义时果断崩溃
        return sqrt(x);
    }
    return sqrt(fabs(x));
}

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

static inline double dd_legendre_1d(int i, double x) {
    if (i == 0 || i == 1) return 0.0;
    if (i == 2) return 3.0;
    return 0.0;
}

void precompute_matrices(void) {
    for (int k = 0; k < N_MODE; k++) {
        M_diag_inv[k] = (2.0 * k + 1.0) / 2.0; 
        for (int q = 0; q < N_QUAD; q++) {
            phi_vol[k][q]  = legendre_1d(k, nodes1D[q]);
            dphi_vol[k][q] = d_legendre_1d(k, nodes1D[q]);
        }
        phi_face_L[k] = legendre_1d(k, -1.0);
        phi_face_R[k] = legendre_1d(k, 1.0);
    }
}

/* ============================================================
 * 物理量与通量计算
 * ============================================================ */
static inline double calc_pressure(double rho, double rhou, double E) {
    return (gamma - 1.0) * (E - 0.5 * rhou * rhou / rho);
}

static inline void euler_flux(double U[NUM_VARS], double F[NUM_VARS]) {
    double rho = U[0];
    double rhou = U[1], E = U[2];
    double u = rhou / rho;
    double p = calc_pressure(rho, rhou, E);
    F[0] = rhou;
    F[1] = rhou * u + p;
    F[2] = u * (E + p);
}

static inline double max_speed() {
    double max_v = 0.0;
    for (int i = 0; i < Nx; i++) {
        double rho = Mesh[i].U[0][0];
        double rhou = Mesh[i].U[1][0];
        double E = Mesh[i].U[2][0];
        double u = rhou / rho;
        double p = calc_pressure(rho, rhou, E);
        double c = safe_sqrt(gamma * p / rho);
        max_v = fmax(max_v, fabs(u) + c);
    }
    return (max_v < 1e-9) ? 1.0 : max_v;
}

static inline void llf_flux(double UL[NUM_VARS], double UR[NUM_VARS], double flux_res[NUM_VARS]) {
    double FL[NUM_VARS], FR[NUM_VARS];
    euler_flux(UL, FL);
    euler_flux(UR, FR);
    
    double rhoL = UL[0], uL = UL[1]/rhoL;
    double pL = calc_pressure(rhoL, UL[1], UL[2]);
    
    double rhoR = UR[0], uR = UR[1]/rhoR;
    double pR = calc_pressure(rhoR, UR[1], UR[2]);
    
    double cL = safe_sqrt(gamma * pL / rhoL);
    double cR = safe_sqrt(gamma * pR / rhoR);
    double alpha = fmax(fabs(uL) + cL, fabs(uR) + cR);
    
    for(int v = 0; v < NUM_VARS; v++) {
        flux_res[v] = 0.5 * (FL[v] + FR[v] - alpha * (UR[v] - UL[v]));
    }
}

/* ============================================================
 * 初始条件与边界
 * ============================================================ */
void init_condition(void) {
    double rho_L = 1.0, p_L = 1.0, u_L = 0.0;
    double rho_R = 0.125, p_R = 0.1, u_R = 0.0;
    int mid = Nx / 2; // 隔膜位置通常设在网格中点

    for (int i = 0; i < Nx; i++) {
        // 1. 清零所有模态 (保持高阶模态为0，初始化为常数分布)
        for(int v = 0; v < NUM_VARS; v++) {
            for(int k = 0; k < N_MODE; k++) {
                Mesh[i].U[v][k] = 0.0;
            }
        }

        // 2. 根据位置设置初始物理量
        double current_rho, current_u, current_p;
        
        if (i < mid) {
            current_rho = rho_L;
            current_u   = u_L;
            current_p   = p_L;
        } else {
            current_rho = rho_R;
            current_u   = u_R;
            current_p   = p_R;
        }

        // 3. 转换为守恒变量并赋值给第0阶模态 (均值)
        // U[0]: 密度 rho
        Mesh[i].U[0][0] = current_rho;
        
        // U[1]: 动量 m = rho * u
        Mesh[i].U[1][0] = current_rho * current_u;
        
        // U[2]: 总能量 E = p/(gamma-1) + 0.5 * rho * u^2
        Mesh[i].U[2][0] = current_p / (gamma - 1.0) + 0.5 * current_rho * current_u * current_u;
    }
}

void apply_ghost_cells(void) {
    memcpy(Ghost_left.U, Mesh[0].U, NUM_VARS * N_MODE * sizeof(double));
    memcpy(Ghost_right.U, Mesh[Nx - 1].U, NUM_VARS * N_MODE * sizeof(double));
}

void boundary_value(void) {
    for (int idx = 0; idx < Nx; idx++) {
        for (int v = 0; v < NUM_VARS; v++) {
            double vl = 0.0, vr = 0.0;
            for (int k = 0; k < N_MODE; k++) {
                vl += phi_face_L[k] * Mesh[idx].U[v][k];
                vr += phi_face_R[k] * Mesh[idx].U[v][k];
            }
            Mesh[idx].face_L[v] = vl;
            Mesh[idx].face_R[v] = vr;
        }
    }
    for (int v = 0; v < NUM_VARS; v++) {
        double l_vr = 0.0, r_vl = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            l_vr += phi_face_R[k] * Ghost_left.U[v][k];
            r_vl += phi_face_L[k] * Ghost_right.U[v][k];
        }
        Ghost_left.face_R[v] = l_vr;
        Ghost_right.face_L[v] = r_vl;
    }
}

/* ============================================================
 * OFDG 阻尼特征矩阵构造与计算
 * ============================================================ */
static void build_eigen_1d(double rho, double u, double p, double L[NUM_VARS][NUM_VARS]) {
    double c = safe_sqrt(gamma * p / rho);
    if (c < 1e-14 && c >= 0.0) c = 1e-14; // 防止除以真零

    double B1 = gamma - 1.0;
    double inv_c = 1.0 / c;
    double prefix = B1 * inv_c; 

    // 第一行
    L[0][0] = prefix * (0.5 * u * c + 0.25 * B1 * u * u);
    L[0][1] = prefix * (-0.5 * B1 * u - 0.5 * c);
    L[0][2] = prefix * (0.5 * B1);

    // 第二行
    L[1][0] = prefix * (c * c - 0.5 * B1 * u * u);
    L[1][1] = prefix * (B1 * u);
    L[1][2] = prefix * (1.0 - gamma); 

    // 第三行
    L[2][0] = prefix * (-0.5 * u * c + 0.25 * B1 * u * u);
    L[2][1] = prefix * (-0.5 * B1 * u + 0.5 * c);
    L[2][2] = prefix * (0.5 * B1);
}

static void eval_deriv(Element *cell, double xi, double out[3][NUM_VARS]) {
    double scale1 = 2.0 / dx;
    double scale2 = 4.0 / (dx * dx);

    for (int v = 0; v < NUM_VARS; v++) {
        out[0][v] = cell->U[v][0]*legendre_1d(0, xi) + cell->U[v][1]*legendre_1d(1, xi) + cell->U[v][2]*legendre_1d(2, xi);
        out[1][v] = (cell->U[v][1]*d_legendre_1d(1, xi) + cell->U[v][2]*d_legendre_1d(2, xi)) * scale1;
        out[2][v] = (cell->U[v][2]*dd_legendre_1d(2, xi)) * scale2;
    }
}

void compute_damp_coefficients(void) {
    double jump_sq[Nx + 1][3][NUM_VARS]; 
    
    for (int i = 0; i <= Nx; i++) {
        Element *L_cell = (i == 0) ? &Ghost_left : &Mesh[i - 1];
        Element *R_cell = (i == Nx) ? &Ghost_right : &Mesh[i];

        double dL[3][NUM_VARS], dR[3][NUM_VARS];
        eval_deriv(L_cell, 1.0, dL);
        eval_deriv(R_cell, -1.0, dR);

        double rhoL = L_cell->U[0][0];
        double uL   = L_cell->U[1][0] / rhoL;
        double pL   = calc_pressure(rhoL, L_cell->U[1][0], L_cell->U[2][0]);
        double HL   = (L_cell->U[2][0] + pL) / rhoL; 

        double rhoR = R_cell->U[0][0];
        double uR   = R_cell->U[1][0] / rhoR;
        double pR   = calc_pressure(rhoR, R_cell->U[1][0], R_cell->U[2][0]);
        double HR   = (R_cell->U[2][0] + pR) / rhoR; 

        double w = safe_sqrt(rhoR / rhoL);

        double u_hat = (uL + w * uR) / (1.0 + w);
        double H_hat = (HL + w * HR) / (1.0 + w);
        double rho_hat = safe_sqrt(rhoL * rhoR); 

        double B1 = gamma - 1.0;
        double c2_hat = B1 * (H_hat - 0.5 * u_hat * u_hat);
        double p_hat = (rho_hat * c2_hat) / gamma;

        double L_mat[NUM_VARS][NUM_VARS];
        build_eigen_1d(rho_hat, u_hat, p_hat, L_mat);
        
        for (int l = 0; l <= 2; l++) {
            for (int s = 0; s < NUM_VARS; s++) {
                double char_jump = 0.0;
                for (int v = 0; v < NUM_VARS; v++) {
                    char_jump += L_mat[s][v] * (dR[l][v] - dL[l][v]);
                }
                jump_sq[i][l][s] = char_jump * char_jump;
            }
        }
    }
    
    for (int i = 0; i < Nx; i++) {
        for (int l = 0; l <= 2; l++) {
            double max_val = 0.0;
            for (int s = 0; s < NUM_VARS; s++) {
                double val = safe_sqrt(jump_sq[i][l][s] + jump_sq[i+1][l][s]);
                if (val > max_val) max_val = val;
            }
            double fact = (l == 0) ? 1.0 : (l == 1 ? 1.0 : 2.0);
            double h_pow = (l == 0) ? 1.0 : (l == 1 ? dx : dx * dx);
            damp_local[i][l] = (2.0 * (2*l + 1) / 3.0) * (h_pow / fact) * max_val;
        }
    }
}

void apply_ofdg_damping(double RHS[Nx][NUM_VARS][N_MODE]) {
    for (int i = 0; i < Nx; i++) {
        double s1 = (damp_local[i][0] + damp_local[i][1]) / dx;
        double s2 = (damp_local[i][0] + damp_local[i][1] + damp_local[i][2]) / dx;
        
        for (int v = 0; v < NUM_VARS; v++) {
            RHS[i][v][1] -= s1 * Mesh[i].U[v][1];
            RHS[i][v][2] -= s2 * Mesh[i].U[v][2];
        }
    }
}

/* ============================================================
 * DG RHS 
 * ============================================================ */
void compute_rhs(double RHS[Nx][NUM_VARS][N_MODE]) {
    apply_ghost_cells();  
    boundary_value();     

    for (int i = 0; i < Nx; i++) {
        Element *cell = &Mesh[i];
        double Vol_Int[NUM_VARS][N_MODE] = {{0}};
        for (int q = 0; q < N_QUAD; q++) {
            double U_phys[NUM_VARS] = {0};
            for (int k = 0; k < N_MODE; k++)
                for (int v = 0; v < NUM_VARS; v++) U_phys[v] += cell->U[v][k] * phi_vol[k][q];
            
            double F_val[NUM_VARS];
            euler_flux(U_phys, F_val);

            for (int k = 0; k < N_MODE; k++) {
                for (int v = 0; v < NUM_VARS; v++) {
                    Vol_Int[v][k] += weights1D[q] * F_val[v] * dphi_vol[k][q];
                }
            }
        }

        double fL[NUM_VARS], fR[NUM_VARS];
        double U_left_neb[NUM_VARS], U_right_neb[NUM_VARS];
        
        if (i == 0) memcpy(U_left_neb, Ghost_left.face_R, sizeof(U_left_neb));
        else        memcpy(U_left_neb, Mesh[i-1].face_R, sizeof(U_left_neb));
        
        if (i == Nx - 1) memcpy(U_right_neb, Ghost_right.face_L, sizeof(U_right_neb));
        else             memcpy(U_right_neb, Mesh[i+1].face_L, sizeof(U_right_neb));

        llf_flux(U_left_neb, cell->face_L, fL);
        llf_flux(cell->face_R, U_right_neb, fR);

        for (int v = 0; v < NUM_VARS; v++) {
            for (int k = 0; k < N_MODE; k++) {
                double Surf_Int = fR[v] * phi_face_R[k] - fL[v] * phi_face_L[k];
                RHS[i][v][k] = (Vol_Int[v][k] - Surf_Int) * M_diag_inv[k] * (2.0 / dx);
            }
        }
    }
}

/* ============================================================
 * 时间推进 (严密修正式 Modified Exp-RK3)
 * ============================================================ */
static double U0[Nx][NUM_VARS][N_MODE], U1[Nx][NUM_VARS][N_MODE];
static double U2[Nx][NUM_VARS][N_MODE], RHS_buf[Nx][NUM_VARS][N_MODE];

void rk3_step(double dt) {
    for (int i = 0; i < Nx; i++) memcpy(U0[i], Mesh[i].U, sizeof(U0[i]));
    
    compute_damp_coefficients(); 
    
    double a0_global = 0.0;
    for (int i = 0; i < Nx; i++) {
        double val = (damp_local[i][0] + damp_local[i][1] + damp_local[i][2]) / dx;
        if (val > a0_global) {
            a0_global = val;
        }
    }
    
    // 【Stage 1】
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 

    double z1 = a0_global * dt;
    double s1_const = 1.0 + z1 + 0.5*z1*z1 + (1.0/6.0)*z1*z1*z1;

    for (int i = 0; i < Nx; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            U1[i][v][0] = U0[i][v][0] + dt * RHS_buf[i][v][0];
            for (int q = 1; q < N_MODE; q++)
                U1[i][v][q] = (U0[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U0[i][v][q]) / s1_const;
        }
    }
    for (int i = 0; i < Nx; i++) memcpy(Mesh[i].U, U1[i], sizeof(U1[i]));

    // 【Stage 2】
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 
    
    double z2 = a0_global * dt; 
    double s1_stage2 = 1.0 + z2 + 0.5*z2*z2 + (1.0/6.0)*z2*z2*z2;
    double s2_stage2 = 1.0 + 0.5*z2 + 0.125*z2*z2 + (1.0/48.0)*z2*z2*z2;

    for (int i = 0; i < Nx; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            U2[i][v][0] = 0.75 * U0[i][v][0] + 0.25 * (U1[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++)
                U2[i][v][q] = (0.75 / s2_stage2) * U0[i][v][q] + (s1_stage2 / (4.0 * s2_stage2)) * (U1[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U1[i][v][q]);
        }
    }
    for (int i = 0; i < Nx; i++) memcpy(Mesh[i].U, U2[i], sizeof(U2[i]));

    // 【Stage 3】
    compute_rhs(RHS_buf);
    apply_ofdg_damping(RHS_buf); 
    
    for (int i = 0; i < Nx; i++) {
        for (int v = 0; v < NUM_VARS; v++) {
            Mesh[i].U[v][0] = (1.0 / 3.0) * U0[i][v][0] + (2.0 / 3.0) * (U2[i][v][0] + dt * RHS_buf[i][v][0]);
            for (int q = 1; q < N_MODE; q++)
                Mesh[i].U[v][q] = (1.0 / (3.0 * s1_stage2)) * U0[i][v][q] + (2.0 * s2_stage2 / (3.0 * s1_stage2)) * (U2[i][v][q] + dt * RHS_buf[i][v][q] + dt * a0_global * U2[i][v][q]);
        }
    }
}

void output_results(void) {
    FILE *fp = fopen("sod_1d.dat", "w");
    fprintf(fp, "VARIABLES = \"X\", \"Rho\", \"U\", \"P\"\n");
    for (int i = 0; i < Nx; i++) {
        double xc = xL + (i + 0.5) * dx;
        double rho = Mesh[i].U[0][0];
        double u = Mesh[i].U[1][0] / rho;
        double p = calc_pressure(rho, Mesh[i].U[1][0], Mesh[i].U[2][0]);
        fprintf(fp, "%lf %lf %lf %lf\n", xc, rho, u, p);
    }
    fclose(fp);
}

int main(void) {
    precompute_matrices();
    init_condition();
    
    double t = 0.0;
    int nit = 0;
    printf("%-10s  %-14s %-14s\n", "Step", "Time", "dt");
    printf("----------------------------------------\n");

    while (t < T_END) {
        double dt = CFL * dx / max_speed();
        if (t + dt > T_END) dt = T_END - t;

        rk3_step(dt);
        t += dt;
        nit++;
        
        if (nit % 10 == 0 || t >= T_END) {
            printf("%-10d  %-14.6e %-14.6e\n", nit, t, dt);
        }
    }
    output_results();
    return 0;
}