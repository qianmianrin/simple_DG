#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define M_PI 3.14159265358979323846
#define Nx 480
#define Ny 160      // 前台阶问题计算域网格
#define xL 0.0
#define xR 3.0
#define yL 0.0
#define yR 1.0
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma_gas 1.4      

#define NUM_VARS 4 
#define STEP_X 0.6  // 阶梯 x 坐标
#define STEP_Y 0.2  // 阶梯 y 坐标

#define CFL   0.1
#define T_END 4.0          // 前台阶问题终止时间

// ================== P^3 升级核心宏定义 ==================
#define N_MODE 10        
#define N_QUAD 16        
#define N_FACE_PTS 4     
#define N_ZS_PTS 32      
#define N_DERIV 10       

// 核心台阶屏蔽宏
#define IS_IN_STEP(ii, jj) (((ii) >= 0 && (ii) < Nx && (jj) >= 0 && (jj) < Ny) ? ((xL + ((ii) + 0.5) * dx > STEP_X - 1e-6) && (yL + ((jj) + 0.5) * dy < STEP_Y - 1e-6)) : false)

typedef struct {
    double U[NUM_VARS][N_MODE]; 
    double face_top[NUM_VARS][N_FACE_PTS];
    double face_bottom[NUM_VARS][N_FACE_PTS];
    double face_left[NUM_VARS][N_FACE_PTS];
    double face_right[NUM_VARS][N_FACE_PTS];
} Element;

Element Mesh_h[Nx * Ny];

// ================== 积分点与映射定义 (4点精度) ==================
static const double nodes_G_h[N_FACE_PTS]   = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
static const double weights_G_h[N_FACE_PTS] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};

static const double nodes_GL_h[N_FACE_PTS]   = {-1.0, -0.4472135954999579, 0.4472135954999579, 1.0};
static const double weights_GL_h[N_FACE_PTS] = {0.1666666666666666, 0.8333333333333333, 0.8333333333333333, 0.1666666666666666};

static const int mk_map_h[N_MODE] = {0, 1, 0, 2, 1, 0, 3, 2, 1, 0};
static const int nk_map_h[N_MODE] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};

// --- 新增：10阶偏导数对应的 x 与 y 求导次数映射 ---
// 导数顺序：0:val, 1:x, 2:y, 3:xx, 4:xy, 5:yy, 6:xxx, 7:xxy, 8:xyy, 9:yyy
static const int nx_map_deriv_h[N_DERIV] = {0, 1, 0, 2, 1, 0, 3, 2, 1, 0};
static const int ny_map_deriv_h[N_DERIV] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};

double r_quad_h[N_QUAD], s_quad_h[N_QUAD], w_quad_h[N_QUAD]; 
double M_diag_inv_h[N_MODE];          
double phi_vol_h[N_MODE][N_QUAD];      
double dphi_dr_vol_h[N_MODE][N_QUAD];  
double dphi_ds_vol_h[N_MODE][N_QUAD];  
double phi_face_T_h[N_MODE][N_FACE_PTS], phi_face_B_h[N_MODE][N_FACE_PTS];
double phi_face_L_h[N_MODE][N_FACE_PTS], phi_face_R_h[N_MODE][N_FACE_PTS];
double phi_ZS_h[N_MODE][N_ZS_PTS]; 

__constant__ double d_nodes_G[N_FACE_PTS];
__constant__ double d_weights_G[N_FACE_PTS];
__constant__ int d_mk_map[N_MODE];
__constant__ int d_nk_map[N_MODE];
__constant__ int d_nx_map_deriv[N_DERIV];
__constant__ int d_ny_map_deriv[N_DERIV];

__constant__ double d_r_quad[N_QUAD], d_s_quad[N_QUAD], d_w_quad[N_QUAD];
__constant__ double d_M_diag_inv[N_MODE];
__constant__ double d_phi_vol[N_MODE][N_QUAD];
__constant__ double d_dphi_dr_vol[N_MODE][N_QUAD];
__constant__ double d_dphi_ds_vol[N_MODE][N_QUAD];
__constant__ double d_phi_face_T[N_MODE][N_FACE_PTS], d_phi_face_B[N_MODE][N_FACE_PTS];
__constant__ double d_phi_face_L[N_MODE][N_FACE_PTS], d_phi_face_R[N_MODE][N_FACE_PTS];
__constant__ double d_phi_ZS[N_MODE][N_ZS_PTS];

// ================== 物理逻辑与基础函数 ==================
__device__ static inline double safe_sqrt(double x) {
  return sqrt(x);
}

__device__ static inline double calc_pressure(double rho, double rhou, double rhov, double E) {
    return (gamma_gas - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}

__device__ static inline void euler_flux(double U[NUM_VARS], double F[NUM_VARS], double G[NUM_VARS]) {
    double rho = U[0], rhou = U[1], rhov = U[2], E = U[3];
    double u = rhou / rho, v = rhov / rho;
    double p = calc_pressure(rho, rhou, rhov, E);

    F[0] = rhou;       F[1] = rhou * u + p; F[2] = rhou * v;       F[3] = u * (E + p);
    G[0] = rhov;       G[1] = rhou * v;     G[2] = rhov * v + p;   G[3] = v * (E + p);
}

__device__ static inline void llf_flux_vector(double UL[NUM_VARS], double UR[NUM_VARS], double nx, double ny, double flux_res[NUM_VARS]) {
    double FL[NUM_VARS], GL[NUM_VARS], FR[NUM_VARS], GR[NUM_VARS];
    euler_flux(UL, FL, GL);
    euler_flux(UR, FR, GR);
    
    double rhoL = UL[0]; double uL = UL[1]/rhoL, vL = UL[2]/rhoL;
    double pL = calc_pressure(rhoL, UL[1], UL[2], UL[3]);
    double cL = safe_sqrt(gamma_gas * pL / rhoL);
    double unL = uL * nx + vL * ny;
    
    double rhoR = UR[0]; double uR = UR[1]/rhoR, vR = UR[2]/rhoR;
    double pR = calc_pressure(rhoR, UR[1], UR[2], UR[3]);
    double cR = safe_sqrt(gamma_gas * pR / rhoR);
    double unR = uR * nx + vR * ny;
    
    double alpha = fmax(fabs(unL) + cL, fabs(unR) + cR);
    
    for(int v = 0; v < NUM_VARS; v++) {
        double flux_n_L = FL[v] * nx + GL[v] * ny;
        double flux_n_R = FR[v] * nx + GR[v] * ny;
        flux_res[v] = 0.5 * (flux_n_L + flux_n_R - alpha * (UR[v] - UL[v]));
    }
}

__device__ static void eval_legendre_basis_1d(double x, double P[4], double dP[4], double ddP[4], double dddP[4]) {
    P[0] = 1.0; dP[0] = 0.0; ddP[0] = 0.0; dddP[0] = 0.0;
    P[1] = x;   dP[1] = 1.0; ddP[1] = 0.0; dddP[1] = 0.0;
    P[2] = 0.5 * (3.0 * x * x - 1.0); dP[2] = 3.0 * x; ddP[2] = 3.0; dddP[2] = 0.0;
    P[3] = 0.5 * (5.0 * x * x * x - 3.0 * x); dP[3] = 1.5 * (5.0 * x * x - 1.0); ddP[3] = 15.0 * x; dddP[3] = 15.0;
}

__device__ static void eval_element_derivatives(Element *cell, double xi, double eta, double deriv_out[N_DERIV][NUM_VARS]) {
    double P_xi[4], dP_xi[4], ddP_xi[4], dddP_xi[4];
    double P_eta[4], dP_eta[4], ddP_eta[4], dddP_eta[4];

    eval_legendre_basis_1d(xi, P_xi, dP_xi, ddP_xi, dddP_xi);
    eval_legendre_basis_1d(eta, P_eta, dP_eta, ddP_eta, dddP_eta);

    double d_dx = 2.0 / dx;
    double d_dy = 2.0 / dy;

    for (int v = 0; v < NUM_VARS; v++) {
        for (int d = 0; d < N_DERIV; d++) deriv_out[d][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            int mk = d_mk_map[k], nk = d_nk_map[k];
            double u_k = cell->U[v][k]; 
            deriv_out[0][v] += u_k * P_xi[mk] * P_eta[nk];
            deriv_out[1][v] += u_k * dP_xi[mk] * P_eta[nk] * d_dx;
            deriv_out[2][v] += u_k * P_xi[mk] * dP_eta[nk] * d_dy;
            deriv_out[3][v] += u_k * ddP_xi[mk] * P_eta[nk] * (d_dx * d_dx);
            deriv_out[4][v] += u_k * dP_xi[mk] * dP_eta[nk] * (d_dx * d_dy);
            deriv_out[5][v] += u_k * P_xi[mk] * ddP_eta[nk] * (d_dy * d_dy);
            deriv_out[6][v] += u_k * dddP_xi[mk] * P_eta[nk] * (d_dx * d_dx * d_dx);
            deriv_out[7][v] += u_k * ddP_xi[mk] * dP_eta[nk] * (d_dx * d_dx * d_dy);
            deriv_out[8][v] += u_k * dP_xi[mk] * ddP_eta[nk] * (d_dx * d_dy * d_dy);
            deriv_out[9][v] += u_k * P_xi[mk] * dddP_eta[nk] * (d_dy * d_dy * d_dy);
        }
    }
}

__device__ void compute_face_values_device(Element *cell) {
    for (int var = 0; var < NUM_VARS; var++) {
        for (int p = 0; p < N_FACE_PTS; p++) {
            double vt = 0.0, vb = 0.0, vl = 0.0, vr = 0.0;
            for (int k = 0; k < N_MODE; k++) {
                double mode_val = cell->U[var][k];
                vt += d_phi_face_T[k][p] * mode_val;
                vb += d_phi_face_B[k][p] * mode_val;
                vl += d_phi_face_L[k][p] * mode_val;
                vr += d_phi_face_R[k][p] * mode_val;
            }
            cell->face_top[var][p]    = vt;
            cell->face_bottom[var][p] = vb;
            cell->face_left[var][p]   = vl;
            cell->face_right[var][p]  = vr;
        }
    }
}

__global__ void apply_ghost_cells_kernel(Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                         Element *d_Ghost_left, Element *d_Ghost_right, double current_time) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Inflow State (gamma, 3, 0, 1)
    const double Q_in[4] = {gamma_gas, gamma_gas * 3.0, 0.0, 1.0/(gamma_gas-1.0) + 0.5 * gamma_gas * 9.0};

    if (id < Nx) {
        // --- 底部边界 (y=0) --- 统一作为反射壁面处理
        for (int k = 0; k < N_MODE; k++) {
            int nk = d_nk_map[k]; 
            double sign = (nk % 2 == 0) ? 1.0 : -1.0; 
            d_Ghost_bottom[id].U[0][k] =  sign * d_Mesh[0 * Nx + id].U[0][k];
            d_Ghost_bottom[id].U[1][k] =  sign * d_Mesh[0 * Nx + id].U[1][k];
            d_Ghost_bottom[id].U[2][k] = -sign * d_Mesh[0 * Nx + id].U[2][k]; // 反转 V
            d_Ghost_bottom[id].U[3][k] =  sign * d_Mesh[0 * Nx + id].U[3][k];
        }
        compute_face_values_device(&d_Ghost_bottom[id]);

        // --- 顶部边界 (y=1) --- 统一作为反射壁面处理
        for (int k = 0; k < N_MODE; k++) {
            int nk = d_nk_map[k]; 
            double sign = (nk % 2 == 0) ? 1.0 : -1.0; 
            d_Ghost_top[id].U[0][k] =  sign * d_Mesh[(Ny - 1) * Nx + id].U[0][k];
            d_Ghost_top[id].U[1][k] =  sign * d_Mesh[(Ny - 1) * Nx + id].U[1][k];
            d_Ghost_top[id].U[2][k] = -sign * d_Mesh[(Ny - 1) * Nx + id].U[2][k]; // 反转 V
            d_Ghost_top[id].U[3][k] =  sign * d_Mesh[(Ny - 1) * Nx + id].U[3][k];
        }
        compute_face_values_device(&d_Ghost_top[id]);
    }
    
    if (id < Ny) {
        // --- 左边界 (x=0) --- Inflow
        for (int v = 0; v < NUM_VARS; v++) {
            d_Ghost_left[id].U[v][0] = Q_in[v];
            for (int k = 1; k < N_MODE; k++) d_Ghost_left[id].U[v][k] = 0.0;
        }
        compute_face_values_device(&d_Ghost_left[id]);

        // --- 右边界 (x=3) --- Outflow (零梯度外推)
        for (int v = 0; v < NUM_VARS; v++) {
            for (int k = 0; k < N_MODE; k++) {
                d_Ghost_right[id].U[v][k] = d_Mesh[id * Nx + (Nx - 1)].U[v][k];
            }
        }
        compute_face_values_device(&d_Ghost_right[id]);
    }
}

__global__ void compute_boundary_faces_kernel(Element *d_Mesh) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii < Nx && jj < Ny) {
        if (IS_IN_STEP(ii, jj)) return; // 屏蔽
        compute_face_values_device(&d_Mesh[jj * Nx + ii]);
    }
}

__device__ void get_vertex_deriv_device(int ii, int jj, int vert, double out[N_DERIV][NUM_VARS], 
                                        Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                        Element *d_Ghost_left, Element *d_Ghost_right,
                                        double d_cell_vertex_derivs[Nx*Ny][4][N_DERIV][NUM_VARS]) {
    // 若请求的是台阶内部坐标，直接返回 0（外部会做精确镜像拦截）
    if (IS_IN_STEP(ii, jj)) {
        for(int d=0; d<N_DERIV; d++) for(int v=0; v<NUM_VARS; v++) out[d][v] = 0.0;
        return;
    }

    if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) {
        for(int d=0; d<N_DERIV; d++) for(int v=0; v<NUM_VARS; v++)
            out[d][v] = d_cell_vertex_derivs[jj * Nx + ii][vert][d][v];
    } else {
        Element *ghost = NULL;
        if (jj < 0) ghost = &d_Ghost_bottom[ii];
        else if (jj >= Ny) ghost = &d_Ghost_top[ii];
        else if (ii < 0) ghost = &d_Ghost_left[jj];
        else if (ii >= Nx) ghost = &d_Ghost_right[jj];
        
        double xi = (vert == 1 || vert == 2) ? 1.0 : -1.0;
        double eta = (vert == 2 || vert == 3) ? 1.0 : -1.0;
        eval_element_derivatives(ghost, xi, eta, out);
    }
}

__global__ void precompute_vertex_derivs_kernel(Element *d_Mesh, double d_cell_vertex_derivs[Nx*Ny][4][N_DERIV][NUM_VARS]) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii < Nx && jj < Ny) {
        if (IS_IN_STEP(ii, jj)) return; // 屏蔽
        int idx = jj * Nx + ii;
        eval_element_derivatives(&d_Mesh[idx], -1.0, -1.0, d_cell_vertex_derivs[idx][0]);
        eval_element_derivatives(&d_Mesh[idx],  1.0, -1.0, d_cell_vertex_derivs[idx][1]);
        eval_element_derivatives(&d_Mesh[idx],  1.0,  1.0, d_cell_vertex_derivs[idx][2]);
        eval_element_derivatives(&d_Mesh[idx], -1.0,  1.0, d_cell_vertex_derivs[idx][3]);
    }
}

__global__ void compute_damp_coeffs_kernel(Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                           Element *d_Ghost_left, Element *d_Ghost_right,
                                           double d_cell_vertex_derivs[Nx*Ny][4][N_DERIV][NUM_VARS],
                                           double d_damp_local[Nx*Ny][4]) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= Nx || jj >= Ny) return;
    if (IS_IN_STEP(ii, jj)) return; // 屏蔽台阶内部
    int idx = jj * Nx + ii;

    Element *C_cell = &d_Mesh[idx];
    double deriv_curr[4][N_DERIV][NUM_VARS], deriv_B[2][N_DERIV][NUM_VARS], deriv_T[2][N_DERIV][NUM_VARS], deriv_L[2][N_DERIV][NUM_VARS], deriv_R[2][N_DERIV][NUM_VARS];    

    get_vertex_deriv_device(ii, jj, 0, deriv_curr[0], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj, 1, deriv_curr[1], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj, 2, deriv_curr[2], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj, 3, deriv_curr[3], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 

    get_vertex_deriv_device(ii, jj - 1, 3, deriv_B[0], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj - 1, 2, deriv_B[1], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj + 1, 0, deriv_T[0], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii, jj + 1, 1, deriv_T[1], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii - 1, jj, 1, deriv_L[0], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii - 1, jj, 2, deriv_L[1], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii + 1, jj, 0, deriv_R[0], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 
    get_vertex_deriv_device(ii + 1, jj, 3, deriv_R[1], d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs); 

    // 【核心修复】：对于台阶造成的内部固壁边界，动态构建精确的数学镜像导数
    bool L_is_step = IS_IN_STEP(ii-1, jj);
    bool R_is_step = IS_IN_STEP(ii+1, jj);
    bool B_is_step = IS_IN_STEP(ii, jj-1); 
    bool T_is_step = IS_IN_STEP(ii, jj+1);

    double jump[4][2][N_DERIV][NUM_VARS];
    for(int d=0; d<N_DERIV; d++) {
        // 根据导数中包含的 x 与 y 的求导阶次计算奇偶翻转符号
        double sign_nx = (d_nx_map_deriv[d] % 2 == 0) ? 1.0 : -1.0;
        double sign_ny = (d_ny_map_deriv[d] % 2 == 0) ? 1.0 : -1.0;

        for(int v=0; v<NUM_VARS; v++) {
            double v_flip_x = (v == 1) ? -1.0 : 1.0; // 垂直壁面反转 U 动量
            double v_flip_y = (v == 2) ? -1.0 : 1.0; // 水平壁面反转 V 动量
            
            // 构建台阶邻居虚拟导数
            double dL0 = L_is_step ? (deriv_curr[0][d][v] * sign_nx * v_flip_x) : deriv_L[0][d][v];
            double dL1 = L_is_step ? (deriv_curr[3][d][v] * sign_nx * v_flip_x) : deriv_L[1][d][v];
            double dR0 = R_is_step ? (deriv_curr[1][d][v] * sign_nx * v_flip_x) : deriv_R[0][d][v];
            double dR1 = R_is_step ? (deriv_curr[2][d][v] * sign_nx * v_flip_x) : deriv_R[1][d][v];
            
            double dB0 = B_is_step ? (deriv_curr[0][d][v] * sign_ny * v_flip_y) : deriv_B[0][d][v];
            double dB1 = B_is_step ? (deriv_curr[1][d][v] * sign_ny * v_flip_y) : deriv_B[1][d][v];
            double dT1 = T_is_step ? (deriv_curr[2][d][v] * sign_ny * v_flip_y) : deriv_T[1][d][v];
            double dT0 = T_is_step ? (deriv_curr[3][d][v] * sign_ny * v_flip_y) : deriv_T[0][d][v];

            jump[0][0][d][v] = fabs(deriv_curr[0][d][v] - dL0); 
            jump[0][1][d][v] = fabs(deriv_curr[0][d][v] - dB0);
            jump[1][0][d][v] = fabs(deriv_curr[1][d][v] - dR0); 
            jump[1][1][d][v] = fabs(deriv_curr[1][d][v] - dB1);
            jump[2][0][d][v] = fabs(deriv_curr[2][d][v] - dR1); 
            jump[2][1][d][v] = fabs(deriv_curr[2][d][v] - dT1);
            jump[3][0][d][v] = fabs(deriv_curr[3][d][v] - dL1); 
            jump[3][1][d][v] = fabs(deriv_curr[3][d][v] - dT0);
        }
    }

    double U_mean[NUM_VARS];
    for (int v = 0; v < NUM_VARS; v++) U_mean[v] = C_cell->U[v][0];
    double rho = U_mean[0];
    double rhou = U_mean[1], rhov = U_mean[2], E = U_mean[3];
    double p = calc_pressure(rho, rhou, rhov, E);
    
    double c = safe_sqrt(gamma_gas * p / rho);
    double H = (E + p) / rho;
    double beta_x = fabs(rhou / rho) + c;
    double beta_y = fabs(rhov / rho) + c;

    int order_map[N_DERIV] = {0, 1, 1, 2, 2, 2, 3, 3, 3, 3}; 
    double sum_X[4][NUM_VARS] = {0};
    double sum_Y[4][NUM_VARS] = {0};

    for(int d=0; d<N_DERIV; d++) {
        int m = order_map[d];
        for(int s=0; s<NUM_VARS; s++) {
            sum_X[m][s] += 0.5*(jump[0][0][d][s] + jump[3][0][d][s] + jump[1][0][d][s] + jump[2][0][d][s]);
            sum_Y[m][s] += 0.5*(jump[0][1][d][s] + jump[1][1][d][s] + jump[3][1][d][s] + jump[2][1][d][s]);
        }
    }

    d_damp_local[idx][0] = 0.0;
    for(int l=1; l<=3; l++) {
        double max_val = 0.0;
        for(int s=0; s<NUM_VARS; s++) {
            double S_X = dx * sum_X[0][s];
            double S_Y = dy * sum_Y[0][s];
            if (l >= 1) { S_X += 2.0 * dx * dx * sum_X[1][s]; S_Y += 2.0 * dy * dy * sum_Y[1][s]; }
            if (l >= 2) { S_X += 6.0 * dx * dx * dx * sum_X[2][s]; S_Y += 6.0 * dy * dy * dy * sum_Y[2][s]; }
            if (l >= 3) { S_X += 12.0 * pow(dx, 4) * sum_X[3][s]; S_Y += 12.0 * pow(dy, 4) * sum_Y[3][s]; }
            double val = 0.0;
            if (H > 1e-13) val = (1.0 / H) * ( (beta_x / dx) * S_X + (beta_y / dy) * S_Y );
            if (val > max_val) max_val = val;
        }
        d_damp_local[idx][l] = max_val;
    }
}

__device__ void get_neighbor_face_device(int ii, int jj, int face, double u_plus[NUM_VARS][N_FACE_PTS],
                                         Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                         Element *d_Ghost_left, Element *d_Ghost_right) {
    // 拦截台阶内部，强制供给精确面通量镜像
    bool is_step_top = (face == 0 && IS_IN_STEP(ii, jj-1));
    bool is_step_front = (face == 3 && IS_IN_STEP(ii+1, jj));

    if (is_step_top) {
        for (int k = 0; k < N_FACE_PTS; k++) {
            u_plus[0][k] = d_Mesh[jj * Nx + ii].face_bottom[0][k];
            u_plus[1][k] = d_Mesh[jj * Nx + ii].face_bottom[1][k];
            u_plus[2][k] = -d_Mesh[jj * Nx + ii].face_bottom[2][k]; 
            u_plus[3][k] = d_Mesh[jj * Nx + ii].face_bottom[3][k];
        }
        return;
    }

    if (is_step_front) {
        for (int k = 0; k < N_FACE_PTS; k++) {
            u_plus[0][k] = d_Mesh[jj * Nx + ii].face_right[0][k];
            u_plus[1][k] = -d_Mesh[jj * Nx + ii].face_right[1][k]; 
            u_plus[2][k] = d_Mesh[jj * Nx + ii].face_right[2][k];
            u_plus[3][k] = d_Mesh[jj * Nx + ii].face_right[3][k];
        }
        return;
    }

    // ============ 标准的全局外边界 ============
    for(int var = 0; var < NUM_VARS; var++) {
        if (face == 0) { 
            if (jj == 0) for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Ghost_bottom[ii].face_top[var][k];
            else for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Mesh[(jj - 1) * Nx + ii].face_top[var][k];
        } else if (face == 1) { 
            if (jj == Ny - 1) for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Ghost_top[ii].face_bottom[var][k];
            else for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Mesh[(jj + 1) * Nx + ii].face_bottom[var][k];
        } else if (face == 2) { 
            if (ii == 0) for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Ghost_left[jj].face_right[var][k];
            else for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Mesh[jj * Nx + (ii - 1)].face_right[var][k];
        } else { 
            if (ii == Nx - 1) for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Ghost_right[jj].face_left[var][k];
            else for (int k = 0; k < N_FACE_PTS; k++) u_plus[var][k] = d_Mesh[jj * Nx + (ii + 1)].face_left[var][k];
        }
    }
}

__global__ void compute_rhs_kernel(Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                   Element *d_Ghost_left, Element *d_Ghost_right,
                                   double d_RHS[Nx*Ny][NUM_VARS][N_MODE]) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= Nx || jj >= Ny) return;
    int idx = jj * Nx + ii;

    // === 内部屏蔽 ===
    if (IS_IN_STEP(ii, jj)) {
        for(int v=0; v<NUM_VARS; v++)
            for(int k=0; k<N_MODE; k++) d_RHS[idx][v][k] = 0.0;
        return;
    }

    Element *cell = &d_Mesh[idx];
    double u_plus[NUM_VARS][N_FACE_PTS];
    double Vol_Int[NUM_VARS][N_MODE] = {{0}};
    double Surf_Int[NUM_VARS][N_MODE] = {{0}};

    for (int q = 0; q < N_QUAD; q++) {
        double U_phys[NUM_VARS] = {0};
        for (int k = 0; k < N_MODE; k++)
            for(int v=0; v<NUM_VARS; v++) U_phys[v] += cell->U[v][k] * d_phi_vol[k][q];
        
        double F_val[NUM_VARS], G_val[NUM_VARS];
        euler_flux(U_phys, F_val, G_val);

        for (int k = 0; k < N_MODE; k++) {
            for(int v=0; v<NUM_VARS; v++) {
                Vol_Int[v][k] += d_w_quad[q] * ( F_val[v] * d_dphi_dr_vol[k][q] * (dy / 2.0) + 
                                                 G_val[v] * d_dphi_ds_vol[k][q] * (dx / 2.0) );
            }
        }
    }

    double U_minus[NUM_VARS], U_p[NUM_VARS], num_f[NUM_VARS];

    get_neighbor_face_device(ii, jj, 0, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < N_FACE_PTS; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_bottom[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 0, -1, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_B[k][p] * (dx / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 1, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < N_FACE_PTS; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 0, 1, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_T[k][p] * (dx / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 2, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < N_FACE_PTS; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, -1, 0, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_L[k][p] * (dy / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 3, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < N_FACE_PTS; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 1, 0, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights_G[p] * num_f[v] * d_phi_face_R[k][p] * (dy / 2.0);
    }

    double J = (dx * dy) / 4.0; 
    for(int v=0; v<NUM_VARS; v++) {
        for (int k = 0; k < N_MODE; k++) {
            d_RHS[idx][v][k] = (Vol_Int[v][k] - Surf_Int[v][k]) * d_M_diag_inv[k] / J;
        }
    }
}

__global__ void update_rk3_stage_kernel(Element *d_Mesh, double d_U0[Nx*Ny][NUM_VARS][N_MODE], 
                                        double d_U_prev[Nx*Ny][NUM_VARS][N_MODE], 
                                        double d_RHS[Nx*Ny][NUM_VARS][N_MODE], int stage, double dt) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;
    int ii = id % Nx; int jj = id / Nx;
    if (IS_IN_STEP(ii, jj)) return; // 屏蔽
    for (int v = 0; v < NUM_VARS; v++) {
        for (int k = 0; k < N_MODE; k++) {
            if (stage == 1) d_Mesh[id].U[v][k] = d_U0[id][v][k] + dt * d_RHS[id][v][k];
            else if (stage == 2) d_Mesh[id].U[v][k] = 0.75 * d_U0[id][v][k] + 0.25 * (d_U_prev[id][v][k] + dt * d_RHS[id][v][k]);
            else if (stage == 3) d_Mesh[id].U[v][k] = (1.0 / 3.0) * d_U0[id][v][k] + (2.0 / 3.0) * (d_U_prev[id][v][k] + dt * d_RHS[id][v][k]);
        }
    }
}

__global__ void apply_jump_filter_kernel(Element *d_Mesh, double d_damp_local[Nx*Ny][4], double dt) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;
    int ii = id % Nx; int jj = id / Nx;
    if (IS_IN_STEP(ii, jj)) return; // 屏蔽
    for (int k = 1; k < N_MODE; k++) {
        int l = d_mk_map[k] + d_nk_map[k];
        if (l > 0 && l <= 3) {
            double sigma = d_damp_local[id][l];
            for (int v = 0; v < NUM_VARS; v++) d_Mesh[id].U[v][k] *= exp(-sigma * dt);
        }
    }
}

__global__ void apply_zhang_shu_limiter_kernel(Element *d_Mesh) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;
    int ii = id % Nx; int jj = id / Nx;
    if (IS_IN_STEP(ii, jj)) return; // 屏蔽，防止产生 NaN

    Element *cell = &d_Mesh[id];
    
    double Ubar[NUM_VARS];
    for (int v = 0; v < NUM_VARS; v++) Ubar[v] = cell->U[v][0];

    double p_bar = calc_pressure(Ubar[0], Ubar[1], Ubar[2], Ubar[3]);
    
    const double eps_rho = 1e-13; 
    double eps_p_target = fmax(1e-13, 1e-10 * Ubar[3]);
    double eps_p = fmin(eps_p_target, 0.99 * p_bar);

    double U_test[N_ZS_PTS][NUM_VARS];
    for (int pt = 0; pt < N_ZS_PTS; pt++) {
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            for (int v = 0; v < NUM_VARS; v++) {
                U_test[pt][v] += cell->U[v][k] * d_phi_ZS[k][pt];
            }
        }
    }

    double rho_min = U_test[0][0];
    for (int i = 1; i < N_ZS_PTS; i++) if (U_test[i][0] < rho_min) rho_min = U_test[i][0];

    double theta1 = 1.0;
    if (rho_min < eps_rho) {
        theta1 = (Ubar[0] - eps_rho) / (Ubar[0] - rho_min);
        if (theta1 < 0.0) theta1 = 0.0; if (theta1 > 1.0) theta1 = 1.0;
        for (int v = 0; v < NUM_VARS; v++) for (int k = 1; k < N_MODE; k++) cell->U[v][k] *= theta1;
        for (int i = 0; i < N_ZS_PTS; i++) for (int v = 0; v < NUM_VARS; v++) U_test[i][v] = Ubar[v] + theta1 * (U_test[i][v] - Ubar[v]);
    }

    double theta2 = 1.0;
    for (int i = 0; i < N_ZS_PTS; i++) {
        double p_i = calc_pressure(U_test[i][0], U_test[i][1], U_test[i][2], U_test[i][3]);
        if (p_i < eps_p) {
            double t_L = 0.0, t_R = 1.0;
            for (int iter = 0; iter < 50; iter++) {
                double t_mid = 0.5 * (t_L + t_R);
                double Ut[NUM_VARS];
                for (int v = 0; v < NUM_VARS; v++) Ut[v] = Ubar[v] + t_mid * (U_test[i][v] - Ubar[v]);
                double p_mid = calc_pressure(Ut[0], Ut[1], Ut[2], Ut[3]);
                
                if (p_mid < eps_p) t_R = t_mid; else t_L = t_mid;
            }
            if (t_L < theta2) theta2 = t_L;
        }
    }
    
    if (theta2 < 1.0)
        for (int v = 0; v < NUM_VARS; v++) for (int k = 1; k < N_MODE; k++) cell->U[v][k] *= theta2;
}

__global__ void copy_mesh_to_tmp(Element *d_Mesh, double d_U_tmp[Nx*Ny][NUM_VARS][N_MODE]) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < Nx * Ny) {
        for(int v=0; v<NUM_VARS; v++) for(int k=0; k<N_MODE; k++) d_U_tmp[id][v][k] = d_Mesh[id].U[v][k];
    }
}

__global__ void copy_mesh_to_U0(Element *d_Mesh, double d_U0[Nx*Ny][NUM_VARS][N_MODE]) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < Nx * Ny) {
        for(int v=0; v<NUM_VARS; v++) for(int k=0; k<N_MODE; k++) d_U0[id][v][k] = d_Mesh[id].U[v][k];
    }
}

struct WaveSpeedFunctor {
    Element* mesh;
    WaveSpeedFunctor(Element* _mesh) : mesh(_mesh) {}
    __device__ double operator()(const int& idx) const {
        int ii = idx % Nx;
        int jj = idx / Nx;
        if (IS_IN_STEP(ii, jj)) return 0.0; // 屏蔽

        double rho =mesh[idx].U[0][0]; 
        double rhou = mesh[idx].U[1][0], rhov = mesh[idx].U[2][0], E = mesh[idx].U[3][0];
        double u = rhou / rho, v = rhov / rho;
        double p = calc_pressure(rho, rhou, rhov, E);
        double c = safe_sqrt(gamma_gas * p / rho);
        return fmax(fmax(fabs(u) + c, 0.0), fabs(v) + c);
    }
};

static inline double legendre_1d(int i, double x) {
    if (i == 0) return 1.0; 
    if (i == 1) return x; 
    if (i == 2) return 0.5 * (3.0 * x * x - 1.0); 
    if (i == 3) return 0.5 * (5.0 * x * x * x - 3.0 * x);
    return 0.0;
}

static inline double d_legendre_1d(int i, double x) {
    if (i == 0) return 0.0; 
    if (i == 1) return 1.0; 
    if (i == 2) return 3.0 * x; 
    if (i == 3) return 1.5 * (5.0 * x * x - 1.0);
    return 0.0;
}

void init_quadrature_and_matrices(void) {
    for (int j = 0; j < N_FACE_PTS; j++) for (int i = 0; i < N_FACE_PTS; i++) {
        int idx = j * N_FACE_PTS + i;
        r_quad_h[idx] = nodes_G_h[i]; s_quad_h[idx] = nodes_G_h[j]; w_quad_h[idx] = weights_G_h[i] * weights_G_h[j];
    }
    
    for (int k = 0; k < N_MODE; k++) {
        int mk = mk_map_h[k], nk = nk_map_h[k]; 
        M_diag_inv_h[k] = 1.0 / ((2.0 / (2.0 * mk + 1.0)) * (2.0 / (2.0 * nk + 1.0)));
        for (int q = 0; q < N_QUAD; q++) {
            double r = r_quad_h[q], s = s_quad_h[q];
            phi_vol_h[k][q]     = legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_dr_vol_h[k][q] = d_legendre_1d(mk, r) * legendre_1d(nk, s);
            dphi_ds_vol_h[k][q] = legendre_1d(mk, r) * d_legendre_1d(nk, s);
        }
        for (int p = 0; p < N_FACE_PTS; p++) {
            double np = nodes_G_h[p];
            phi_face_T_h[k][p] = legendre_1d(mk, np) * legendre_1d(nk, 1.0);
            phi_face_B_h[k][p] = legendre_1d(mk, np) * legendre_1d(nk, -1.0);
            phi_face_L_h[k][p] = legendre_1d(mk, -1.0) * legendre_1d(nk, np);
            phi_face_R_h[k][p] = legendre_1d(mk, 1.0) * legendre_1d(nk, np);
        }
    }

    int pt = 0;
    for(int j = 0; j < N_FACE_PTS; j++) { 
        for(int i = 0; i < N_FACE_PTS; i++) {
            double x = nodes_GL_h[i], y = nodes_G_h[j];
            for(int k = 0; k < N_MODE; k++) phi_ZS_h[k][pt] = legendre_1d(mk_map_h[k], x) * legendre_1d(nk_map_h[k], y);
            pt++;
        }
    }
    for(int j = 0; j < N_FACE_PTS; j++) { 
        for(int i = 0; i < N_FACE_PTS; i++) {
            double x = nodes_G_h[i], y = nodes_GL_h[j];
            for(int k = 0; k < N_MODE; k++) phi_ZS_h[k][pt] = legendre_1d(mk_map_h[k], x) * legendre_1d(nk_map_h[k], y);
            pt++;
        }
    }
    for(int j = 0; j < N_FACE_PTS; j++) { 
        for(int i = 0; i < N_FACE_PTS; i++) {
            double x = nodes_G_h[i], y = nodes_G_h[j];
            for(int k = 0; k < N_MODE; k++) {
                phi_ZS_h[k][pt] = legendre_1d(mk_map_h[k], x) * legendre_1d(nk_map_h[k], y);
            }
            pt++; 
        }
    }
}

void init_condition(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh_h[jj * Nx + ii];
            
            // 前台阶问题均匀初始流场
            double rho = gamma_gas;
            double u = 3.0;
            double v = 0.0;
            double p = 1.0;

            double Q[NUM_VARS];
            Q[0] = rho; 
            Q[1] = rho * u; 
            Q[2] = rho * v;
            Q[3] = p / (gamma_gas - 1.0) + 0.5 * rho * (u * u + v * v);

            for(int var = 0; var < NUM_VARS; var++) {
                cell->U[var][0] = Q[var];
                for(int k=1; k<N_MODE; k++) cell->U[var][k] = 0.0;
            }
        }
    }
}

void output_results(double t) {
    FILE *fp = fopen("result2.dat", "w");
    if (fp == NULL) return;
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Rho\", \"U\", \"V\", \"P\"\n");

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            if (IS_IN_STEP(ii, jj)) continue;

            int idx = jj * Nx + ii;
            Element *cell = &Mesh_h[idx]; 
            double xc = xL + (ii + 0.5) * dx;
            double yc = yL + (jj + 0.5) * dy;

            for (int i = 0; i < N_FACE_PTS; i++) {
                for (int j = 0; j < N_FACE_PTS; j++) {
                    double r = nodes_G_h[j], s = nodes_G_h[i]; 
                    double x_phys = xc + (dx / 2.0) * r, y_phys = yc + (dy / 2.0) * s;
                    
                    double U_phys[NUM_VARS] = {0};
                    for (int k = 0; k < N_MODE; k++) {
                        double phi_val = legendre_1d(mk_map_h[k], r) * legendre_1d(nk_map_h[k], s);
                        for (int v = 0; v < NUM_VARS; v++) U_phys[v] += cell->U[v][k] * phi_val;
                    }

                    double rho = U_phys[0];
                    double u = U_phys[1] / rho;
                    double v = U_phys[2] / rho;
                    double p = (gamma_gas - 1.0) * (U_phys[3] - 0.5 * (U_phys[1] * U_phys[1] + U_phys[2] * U_phys[2]) / rho);
                    
                    fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", x_phys, y_phys, rho, u, v, p);
                }
            }
        }
    }
    fclose(fp);
}
__global__ void apply_corner_entropy_fix_kernel(Element *d_Mesh) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id != 0) return; // 仅需单线程处理该奇异点

    // 1. 确定转角位置索引
    int i_corner = round(STEP_X / dx); // x=0.6 对应的单元索引（台阶后的第一个点）
    int j_corner = round(STEP_Y / dy); // y=0.2 对应的单元索引（台阶平面的第一层）

    // 2. 选取参考点 (Reference Zone): 转角左下方的一个单元
    // 对应流体绕过台阶前的上游紧邻单元，未受拐角膨胀扇污染
    int ic_ref = i_corner - 1;
    int jc_ref = j_corner - 1;
    int idx_ref = jc_ref * Nx + ic_ref;

    // 获取参考点的单元平均值 (零阶模态 k=0)
    double rho_ref  = d_Mesh[idx_ref].U[0][0];
    double rhou_ref = d_Mesh[idx_ref].U[1][0];
    double rhov_ref = d_Mesh[idx_ref].U[2][0];
    double E_ref    = d_Mesh[idx_ref].U[3][0];

    double inv_rho_ref = 1.0 / rho_ref;
    double u_ref = rhou_ref * inv_rho_ref;
    double v_ref = rhov_ref * inv_rho_ref;
    double p_ref = (gamma_gas - 1.0) * (E_ref - 0.5 * rho_ref * (u_ref * u_ref + v_ref * v_ref));
    
    // 参考点的单位质量总焓与熵参考值 (P/rho^gamma)
    double H_ref = (E_ref + p_ref) * inv_rho_ref; 
    double entropy_ref = p_ref / pow(rho_ref, gamma_gas); 

    // 3. 修正转角上方的网格
    // 第一层 (j_corner): 修正转角右侧 4 个格点
    // 第二层 (j_corner + 1): 修正转角右侧 2 个格点
    for (int row_off = 0; row_off < 2; row_off++) {
        int j = j_corner + row_off;
        int num_zones = (row_off == 0) ? 4 : 2;

        for (int i_off = 0; i_off < num_zones; i_off++) {
            int i = i_corner + i_off;
            if (i >= Nx || j >= Ny) continue;

            int idx = j * Nx + i;
            
            // 提取目标单元的平均值
            double rho  = d_Mesh[idx].U[0][0];
            double rhou = d_Mesh[idx].U[1][0];
            double rhov = d_Mesh[idx].U[2][0];
            double E    = d_Mesh[idx].U[3][0];

            // A. 计算当前压力 (压力保持不变，因为膨胀波扫过时压力过渡是连续的)
            double inv_rho = 1.0 / rho;
            double u = rhou * inv_rho;
            double v = rhov * inv_rho;
            double p = (gamma_gas - 1.0) * (E - 0.5 * rho * (u * u + v * v));

            // B. 重置密度 (保持与参考点熵相同)
            // rho_new = (p / entropy_ref)^(1/gamma)
            double rho_new = pow(p / entropy_ref, 1.0 / gamma_gas);

            // C. 重置速度模长 (保持总焓相同)
            // H_ref = [gamma*p / ((gamma-1)*rho_new)] + 0.5 * q_new^2
            double enthalpy_new = (gamma_gas * p) / ((gamma_gas - 1.0) * rho_new);
            double q2_new = 2.0 * (H_ref - enthalpy_new);
            if (q2_new < 0.0) q2_new = 0.0; // 能量不足时强制速度为0，防止负开方

            double q_old = sqrt(u * u + v * v);
            if (q_old < 1e-12) q_old = 1e-12; // 防止除以 0
            double scale = sqrt(q2_new) / q_old;

            double u_new = u * scale;
            double v_new = v * scale;

            // D. 写回守恒变量的零阶模态 (均值)
            d_Mesh[idx].U[0][0] = rho_new;
            d_Mesh[idx].U[1][0] = rho_new * u_new;
            d_Mesh[idx].U[2][0] = rho_new * v_new;
            d_Mesh[idx].U[3][0] = p / (gamma_gas - 1.0) + 0.5 * rho_new * (u_new * u_new + v_new * v_new);

            // E. DG 方法的必然步骤：将这 6 个单元的高阶模态抹平
            // 既然其均值已经是通过 0D 热力学状态强行赋予的，旧的高阶梯度必须被抛弃，以阻断数值震荡的繁衍
            for (int k = 1; k < N_MODE; k++) {
                d_Mesh[idx].U[0][k] = 0.0;
                d_Mesh[idx].U[1][k] = 0.0;
                d_Mesh[idx].U[2][k] = 0.0;
                d_Mesh[idx].U[3][k] = 0.0;
            }
        }
    }
}
int main(void) {
    init_quadrature_and_matrices();
    init_condition();

    cudaMemcpyToSymbol(d_nodes_G, nodes_G_h, sizeof(nodes_G_h));
    cudaMemcpyToSymbol(d_weights_G, weights_G_h, sizeof(weights_G_h));
    cudaMemcpyToSymbol(d_mk_map, mk_map_h, sizeof(mk_map_h));
    cudaMemcpyToSymbol(d_nk_map, nk_map_h, sizeof(nk_map_h));
    cudaMemcpyToSymbol(d_nx_map_deriv, nx_map_deriv_h, sizeof(nx_map_deriv_h));
    cudaMemcpyToSymbol(d_ny_map_deriv, ny_map_deriv_h, sizeof(ny_map_deriv_h));
    cudaMemcpyToSymbol(d_r_quad, r_quad_h, sizeof(r_quad_h));
    cudaMemcpyToSymbol(d_s_quad, s_quad_h, sizeof(s_quad_h));
    cudaMemcpyToSymbol(d_w_quad, w_quad_h, sizeof(w_quad_h));
    cudaMemcpyToSymbol(d_M_diag_inv, M_diag_inv_h, sizeof(M_diag_inv_h));
    cudaMemcpyToSymbol(d_phi_vol, phi_vol_h, sizeof(phi_vol_h));
    cudaMemcpyToSymbol(d_dphi_dr_vol, dphi_dr_vol_h, sizeof(dphi_dr_vol_h));
    cudaMemcpyToSymbol(d_dphi_ds_vol, dphi_ds_vol_h, sizeof(dphi_ds_vol_h));
    cudaMemcpyToSymbol(d_phi_face_T, phi_face_T_h, sizeof(phi_face_T_h));
    cudaMemcpyToSymbol(d_phi_face_B, phi_face_B_h, sizeof(phi_face_B_h));
    cudaMemcpyToSymbol(d_phi_face_L, phi_face_L_h, sizeof(phi_face_L_h));
    cudaMemcpyToSymbol(d_phi_face_R, phi_face_R_h, sizeof(phi_face_R_h));
    cudaMemcpyToSymbol(d_phi_ZS, phi_ZS_h, sizeof(phi_ZS_h));

    Element *d_Mesh, *d_Ghost_bottom, *d_Ghost_top, *d_Ghost_left, *d_Ghost_right;
    cudaMalloc(&d_Mesh, Nx * Ny * sizeof(Element));
    cudaMalloc(&d_Ghost_bottom, Nx * sizeof(Element)); cudaMalloc(&d_Ghost_top, Nx * sizeof(Element));
    cudaMalloc(&d_Ghost_left, Ny * sizeof(Element));   cudaMalloc(&d_Ghost_right, Ny * sizeof(Element));

    double (*d_cell_vertex_derivs)[4][N_DERIV][NUM_VARS], (*d_damp_local)[4];
    cudaMalloc(&d_cell_vertex_derivs, Nx * Ny * sizeof(*d_cell_vertex_derivs));
    cudaMalloc(&d_damp_local, Nx * Ny * sizeof(*d_damp_local));

    double (*d_U0)[NUM_VARS][N_MODE], (*d_U_tmp)[NUM_VARS][N_MODE];
    cudaMalloc(&d_U0, Nx * Ny * sizeof(*d_U0)); 
    cudaMalloc(&d_U_tmp, Nx * Ny * sizeof(*d_U_tmp));
    
    double (*d_RHS)[NUM_VARS][N_MODE];
    cudaMalloc(&d_RHS, Nx * Ny * sizeof(*d_RHS));

    cudaMemcpy(d_Mesh, Mesh_h, Nx * Ny * sizeof(Element), cudaMemcpyHostToDevice);

    dim3 blockSize2D(16, 16);
    dim3 gridSize2D((Nx + 15) / 16, (Ny + 15) / 16);
    int blockSize1D = 256;
    int gridSize1D = (Nx * Ny + 255) / 256;
    int gridMaxBoundary = (fmax(Nx, Ny) + blockSize1D - 1) / blockSize1D; 

    thrust::counting_iterator<int> iter(0);

    double current_time = 0.0; 
    int nit = 0;
    printf("%-10s  %-14s %-14s\n", "Step", "Time", "dt");
    printf("------------------------------------------\n");

    while (current_time < T_END) {
        compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
        apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right,current_time);
        
        double max_wave = thrust::transform_reduce(thrust::device, iter, iter + Nx * Ny, WaveSpeedFunctor(d_Mesh), 0.0, thrust::maximum<double>());
        if (max_wave < 1e-9) max_wave = 1.0;
        
        double dt = CFL * fmin(dx, dy) / max_wave;
        if (current_time + dt > T_END) dt = T_END - current_time;

        copy_mesh_to_U0<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0);

      for (int stage = 1; stage <= 3; stage++) {
            if (stage > 1) copy_mesh_to_tmp<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U_tmp);
            
            compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
            apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right,current_time);
            compute_rhs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_RHS);

            if (stage == 1) update_rk3_stage_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0, d_U0, d_RHS, stage, dt);
            else update_rk3_stage_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0, d_U_tmp, d_RHS, stage, dt);

            compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
            apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right,current_time);
            precompute_vertex_derivs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_cell_vertex_derivs);
            
            compute_damp_coeffs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs, d_damp_local);
            
            // 1. 先进行 Jump Filter 阻尼
            apply_jump_filter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_damp_local, dt);

            // ================= 新增：角点物理熵修正 =================
            // 在限制器之前强制剥离数值噪音
            apply_corner_entropy_fix_kernel<<<1, 1>>>(d_Mesh);
            cudaDeviceSynchronize(); 

            // 2. 最后应用保正限制器
            apply_zhang_shu_limiter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh);
        }

        current_time += dt; 
        nit++;
        printf("%-10d  %-14.6e %-14.6e \n", nit, current_time, dt);
    }

    cudaMemcpy(Mesh_h, d_Mesh, Nx * Ny * sizeof(Element), cudaMemcpyDeviceToHost);
    printf("GPU Computation Complete.\n");
    output_results(T_END);

    cudaFree(d_Mesh); cudaFree(d_Ghost_bottom); cudaFree(d_Ghost_top); cudaFree(d_Ghost_left); cudaFree(d_Ghost_right);
    cudaFree(d_cell_vertex_derivs); cudaFree(d_damp_local);
    cudaFree(d_U0); cudaFree(d_U_tmp); cudaFree(d_RHS);
    return 0;
}