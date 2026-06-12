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

/* ============================================================
 * 基本全局参数
 * ============================================================ */

#define Nx 640         // 网格数可自由调整以测试收敛精度
#define Ny 640
#define xL 0.0
#define xR 2.0        // 域宽为 2，正好对应一个完整的正弦周期
#define yL 0.0
#define yR 2.0

#define NUM_VARS 4 
// #define Mach 800.0 // 该宏在此平滑问题中不再需要，可以删除或保留不使用
#define gamma_gas 1.4 // 图中指定 \gamma = 1.4
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)

#define NUM_VARS 4 
#define Mach 800.0 


#define CFL   0.1

#define T_END 1.0    // 图中指定计算到 T = 10


#define N_MODE 6   
#define N_QUAD 9   

typedef struct {
    double U[NUM_VARS][N_MODE]; 
    double face_top[NUM_VARS][3];
    double face_bottom[NUM_VARS][3];
    double face_left[NUM_VARS][3];
    double face_right[NUM_VARS][3];
} Element;

/* ============================================================
 * CPU 端基础数据 (仅用于初始化和输出)
 * ============================================================ */
Element Mesh_h[Nx * Ny];

static const double nodes1D_h[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D_h[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
static const int mk_map_h[N_MODE] = {0, 1, 0, 2, 1, 0};
static const int nk_map_h[N_MODE] = {0, 0, 1, 0, 1, 2};

double r_quad_h[N_QUAD], s_quad_h[N_QUAD], w_quad_h[N_QUAD]; 
double M_diag_inv_h[N_MODE];          
double phi_vol_h[N_MODE][N_QUAD];      
double dphi_dr_vol_h[N_MODE][N_QUAD];  
double dphi_ds_vol_h[N_MODE][N_QUAD];  
double phi_face_T_h[N_MODE][3], phi_face_B_h[N_MODE][3];
double phi_face_L_h[N_MODE][3], phi_face_R_h[N_MODE][3];

/* ============================================================
 * GPU 常量内存 (极速读取)
 * ============================================================ */
__constant__ double d_nodes1D[3];
__constant__ double d_weights1D[3];
__constant__ int d_mk_map[N_MODE];
__constant__ int d_nk_map[N_MODE];

__constant__ double d_r_quad[N_QUAD], d_s_quad[N_QUAD], d_w_quad[N_QUAD];
__constant__ double d_M_diag_inv[N_MODE];
__constant__ double d_phi_vol[N_MODE][N_QUAD];
__constant__ double d_dphi_dr_vol[N_MODE][N_QUAD];
__constant__ double d_dphi_ds_vol[N_MODE][N_QUAD];
__constant__ double d_phi_face_T[N_MODE][3], d_phi_face_B[N_MODE][3];
__constant__ double d_phi_face_L[N_MODE][3], d_phi_face_R[N_MODE][3];

/* ============================================================
 * 物理与数学辅助函数 (Device端)
 * ============================================================ */
__device__ static inline double safe_sqrt(double x) {
    if (x < -1e8) return sqrt(x); 
    return sqrt(fabs(x));
}

__device__ static inline double calc_pressure(double rho, double rhou, double rhov, double E) {
    return (gamma_gas - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}

__device__ static inline double calc_pressure_safe(double rho, double rhou, double rhov, double E) {
    if (rho <= 0.0) return -1.0; 
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
    double FL[NUM_VARS], GL[NUM_VARS];
    double FR[NUM_VARS], GR[NUM_VARS];
    
    euler_flux(UL, FL, GL);
    euler_flux(UR, FR, GR);
    
    double rhoL = UL[0], uL = UL[1]/rhoL, vL = UL[2]/rhoL;
    double pL = calc_pressure(rhoL, UL[1], UL[2], UL[3]);
    double cL = safe_sqrt(gamma_gas * pL / rhoL);
    double unL = uL * nx + vL * ny;
    
    double rhoR = UR[0], uR = UR[1]/rhoR, vR = UR[2]/rhoR;
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


__device__ static void eval_legendre_basis_1d(double x, double P[3], double dP[3], double ddP[3]) {
    P[0] = 1.0; dP[0] = 0.0; ddP[0] = 0.0;
    P[1] = x;   dP[1] = 1.0; ddP[1] = 0.0;
    P[2] = 0.5 * (3.0 * x * x - 1.0); dP[2] = 3.0 * x; ddP[2] = 3.0;
}

__device__ static void eval_element_derivatives(Element *cell, double xi, double eta, double deriv_out[6][NUM_VARS]) {
    double P_xi[3], dP_xi[3], ddP_xi[3];
    double P_eta[3], dP_eta[3], ddP_eta[3];

    eval_legendre_basis_1d(xi, P_xi, dP_xi, ddP_xi);
    eval_legendre_basis_1d(eta, P_eta, dP_eta, ddP_eta);

    double d_dx = 2.0 / dx;
    double d_dy = 2.0 / dy;

    for (int v = 0; v < NUM_VARS; v++) {
        for (int d = 0; d < 6; d++) deriv_out[d][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) {
            int mk = d_mk_map[k], nk = d_nk_map[k];
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


/* ============================================================
 * CUDA Kernels
 * ============================================================ */

__device__ void compute_face_values_device(Element *cell) {
    for (int var = 0; var < NUM_VARS; var++) {
        for (int p = 0; p < 3; p++) {
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
                                         Element *d_Ghost_left, Element *d_Ghost_right) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Y方向边界（上下）的周期性映射
    if (id < Nx) {
        for(int v=0; v<NUM_VARS; v++) {
            for(int k=0; k<N_MODE; k++) {
                // 底部幽灵单元取顶部内部网格的值
                d_Ghost_bottom[id].U[v][k] = d_Mesh[(Ny - 1) * Nx + id].U[v][k];
                // 顶部幽灵单元取底部内部网格的值
                d_Ghost_top[id].U[v][k] = d_Mesh[0 * Nx + id].U[v][k];
            }
        }
        compute_face_values_device(&d_Ghost_bottom[id]);
        compute_face_values_device(&d_Ghost_top[id]);
    }
    
    // X方向边界（左右）的周期性映射
    if (id < Ny) {
        for(int v=0; v<NUM_VARS; v++) {
            for(int k=0; k<N_MODE; k++) {
                // 左侧幽灵单元取右侧内部网格的值
                d_Ghost_left[id].U[v][k] = d_Mesh[id * Nx + (Nx - 1)].U[v][k];
                // 右侧幽灵单元取左侧内部网格的值
                d_Ghost_right[id].U[v][k] = d_Mesh[id * Nx + 0].U[v][k];
            }
        }
        compute_face_values_device(&d_Ghost_left[id]);
        compute_face_values_device(&d_Ghost_right[id]);
    }
}

__global__ void compute_boundary_faces_kernel(Element *d_Mesh) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii < Nx && jj < Ny) {
        compute_face_values_device(&d_Mesh[jj * Nx + ii]);
    }
}

__device__ void get_vertex_deriv_device(int ii, int jj, int vert, double out[6][NUM_VARS], 
                                        Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                        Element *d_Ghost_left, Element *d_Ghost_right,
                                        double d_cell_vertex_derivs[Nx*Ny][4][6][NUM_VARS]) {
    if (ii >= 0 && ii < Nx && jj >= 0 && jj < Ny) {
        for(int d=0; d<6; d++) for(int v=0; v<NUM_VARS; v++)
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

__global__ void precompute_vertex_derivs_kernel(Element *d_Mesh, double d_cell_vertex_derivs[Nx*Ny][4][6][NUM_VARS]) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii < Nx && jj < Ny) {
        int idx = jj * Nx + ii;
        eval_element_derivatives(&d_Mesh[idx], -1.0, -1.0, d_cell_vertex_derivs[idx][0]);
        eval_element_derivatives(&d_Mesh[idx],  1.0, -1.0, d_cell_vertex_derivs[idx][1]);
        eval_element_derivatives(&d_Mesh[idx],  1.0,  1.0, d_cell_vertex_derivs[idx][2]);
        eval_element_derivatives(&d_Mesh[idx], -1.0,  1.0, d_cell_vertex_derivs[idx][3]);
    }
}
__global__ void compute_oedg_damp_kernel(
    Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
    Element *d_Ghost_left, Element *d_Ghost_right,
    double d_cell_vertex_derivs[Nx*Ny][4][6][NUM_VARS],
    double d_damp_local[Nx*Ny][3], 
    double max_dev_rho, double max_dev_rhou, double max_dev_rhov, double max_dev_E) 
{
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= Nx || jj >= Ny) return;
    int idx = jj * Nx + ii;

    Element *C_cell = &d_Mesh[idx];
    
    // 提取顶点导数代码 (同原逻辑保留)...
    double deriv_curr[4][6][NUM_VARS], deriv_B[2][6][NUM_VARS], deriv_T[2][6][NUM_VARS], deriv_L[2][6][NUM_VARS], deriv_R[2][6][NUM_VARS];
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

    double jump[4][2][6][NUM_VARS];
    for(int d=0; d<6; d++) for(int v=0; v<NUM_VARS; v++) {
        jump[0][0][d][v] = deriv_curr[0][d][v] - deriv_L[0][d][v]; jump[0][1][d][v] = deriv_curr[0][d][v] - deriv_B[0][d][v];
        jump[1][0][d][v] = deriv_curr[1][d][v] - deriv_R[0][d][v]; jump[1][1][d][v] = deriv_curr[1][d][v] - deriv_B[1][d][v];
        jump[2][0][d][v] = deriv_curr[2][d][v] - deriv_R[1][d][v]; jump[2][1][d][v] = deriv_curr[2][d][v] - deriv_T[1][d][v];
        jump[3][0][d][v] = deriv_curr[3][d][v] - deriv_L[1][d][v]; jump[3][1][d][v] = deriv_curr[3][d][v] - deriv_T[0][d][v];
    }

    // 局部最大波速 (beta) 计算
    double rho = C_cell->U[0][0], u = C_cell->U[1][0]/rho, v_vel = C_cell->U[2][0]/rho;
    double p = calc_pressure(rho, C_cell->U[1][0], C_cell->U[2][0], C_cell->U[3][0]);
    double c = safe_sqrt(gamma_gas * p / rho);
    double beta_x = fabs(u) + c;
    double beta_y = fabs(v_vel) + c;

    double max_devs[NUM_VARS] = {max_dev_rho, max_dev_rhou, max_dev_rhov, max_dev_E};
    int order_map[6] = {0, 1, 1, 2, 2, 2}; 

    for (int m = 0; m <= 2; m++) {
        if (m == 0) {
            d_damp_local[idx][0] = 0.0; // 数学上严格定义
            continue;
        }
        
        double max_delta_K = 0.0;
        double coef = (2.0 * m + 1.0) / ( (m == 1 ? 6.0 : 12.0)); // k=1或2        
        for (int v_idx = 0; v_idx < NUM_VARS; v_idx++) {
            double denom = fmax(max_devs[v_idx], 1e-12); // 防止全场一致导致的除零
            double sigma_x = 0.0;
            double sigma_y = 0.0;

            for (int d = 0; d < 6; d++) {
                if (order_map[d] != m) continue;
                // 利用梯形法则近似边缘积分 (端点平均)
                double jump_L = 0.5 * (fabs(jump[0][0][d][v_idx]) + fabs(jump[3][0][d][v_idx]));
                double jump_R = 0.5 * (fabs(jump[1][0][d][v_idx]) + fabs(jump[2][0][d][v_idx]));
                double jump_B = 0.5 * (fabs(jump[0][1][d][v_idx]) + fabs(jump[1][1][d][v_idx]));
                double jump_T = 0.5 * (fabs(jump[3][1][d][v_idx]) + fabs(jump[2][1][d][v_idx]));

                sigma_x += (jump_L + jump_R);
                sigma_y += (jump_B + jump_T);
            }

            double h_x_pow = m == 1 ? dx : (dx * dx);
            double h_y_pow = m == 1 ? dy : (dy * dy);
            
            double delta_var = (beta_x * coef * h_x_pow * sigma_x / dx + 
                                beta_y * coef * h_y_pow * sigma_y / dy) / denom;

            if (delta_var > max_delta_K) max_delta_K = delta_var;
        }
        d_damp_local[idx][m] = max_delta_K;
    }
}

__device__ void get_neighbor_face_device(int ii, int jj, int face, double u_plus[NUM_VARS][3],
                                         Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                         Element *d_Ghost_left, Element *d_Ghost_right) {
    for(int var = 0; var < NUM_VARS; var++) {
        if (face == 0) { 
            if (jj == 0) for (int k = 0; k < 3; k++) u_plus[var][k] = d_Ghost_bottom[ii].face_top[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = d_Mesh[(jj - 1) * Nx + ii].face_top[var][k];
        } else if (face == 1) { 
            if (jj == Ny - 1) for (int k = 0; k < 3; k++) u_plus[var][k] = d_Ghost_top[ii].face_bottom[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = d_Mesh[(jj + 1) * Nx + ii].face_bottom[var][k];
        } else if (face == 2) { 
            if (ii == 0) for (int k = 0; k < 3; k++) u_plus[var][k] = d_Ghost_left[jj].face_right[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = d_Mesh[jj * Nx + (ii - 1)].face_right[var][k];
        } else { 
            if (ii == Nx - 1) for (int k = 0; k < 3; k++) u_plus[var][k] = d_Ghost_right[jj].face_left[var][k];
            else for (int k = 0; k < 3; k++) u_plus[var][k] = d_Mesh[jj * Nx + (ii + 1)].face_left[var][k];
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

    Element *cell = &d_Mesh[idx];
    double u_plus[NUM_VARS][3];
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
    for (int p = 0; p < 3; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_bottom[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 0, -1, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights1D[p] * num_f[v] * d_phi_face_B[k][p] * (dx / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 1, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < 3; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 0, 1, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights1D[p] * num_f[v] * d_phi_face_T[k][p] * (dx / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 2, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < 3; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, -1, 0, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights1D[p] * num_f[v] * d_phi_face_L[k][p] * (dy / 2.0);
    }
    
    get_neighbor_face_device(ii, jj, 3, u_plus, d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
    for (int p = 0; p < 3; p++) {
        for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][p]; U_p[v] = u_plus[v][p]; }
        llf_flux_vector(U_minus, U_p, 1, 0, num_f);
        for(int k=0; k<N_MODE; k++) for(int v=0; v<NUM_VARS; v++) 
            Surf_Int[v][k] += d_weights1D[p] * num_f[v] * d_phi_face_R[k][p] * (dy / 2.0);
    }

    double J = (dx * dy) / 4.0; 
    for(int v=0; v<NUM_VARS; v++) {
        for (int k = 0; k < N_MODE; k++) {
            d_RHS[idx][v][k] = (Vol_Int[v][k] - Surf_Int[v][k]) * d_M_diag_inv[k] / J;
        }
    }
}
// 纯物理双曲推进 (标准 SSP-RK3)
__global__ void rk3_physical_stage_kernel(Element *d_Mesh, 
                                          double d_U0[Nx*Ny][NUM_VARS][N_MODE], 
                                          double d_U_prev[Nx*Ny][NUM_VARS][N_MODE], 
                                          double d_RHS[Nx*Ny][NUM_VARS][N_MODE], 
                                          int stage, double dt) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;

    for (int v = 0; v < NUM_VARS; v++) {
        for (int k = 0; k < N_MODE; k++) {
            if (stage == 1) {
                d_Mesh[id].U[v][k] = d_U0[id][v][k] + dt * d_RHS[id][v][k];
            } else if (stage == 2) {
                d_Mesh[id].U[v][k] = 0.75 * d_U0[id][v][k] + 0.25 * (d_U_prev[id][v][k] + dt * d_RHS[id][v][k]);
            } else if (stage == 3) {
                d_Mesh[id].U[v][k] = (1.0 / 3.0) * d_U0[id][v][k] + (2.0 / 3.0) * (d_U_prev[id][v][k] + dt * d_RHS[id][v][k]);
            }
        }
    }
}

// 独立的 OE 模态滤波算子 (指数衰减，绝对无强制截断)
__global__ void apply_oedg_filter_kernel(Element *d_Mesh, double d_damp_local[Nx*Ny][3], double dt) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;

    int order_map[N_MODE] = {0, 1, 1, 2, 2, 2}; 
    
    for (int v = 0; v < NUM_VARS; v++) {
        for (int k = 0; k < N_MODE; k++) {
            int max_m = order_map[k];
            double sum_delta = 0.0;
            for (int m = 0; m <= max_m; m++) {
                sum_delta += d_damp_local[id][m];
            }
            // 利用精确解衰减，m=0 时 sum_delta 为 0，exp(0) = 1，天然守恒。
            d_Mesh[id].U[v][k] *= exp(-dt * sum_delta);
        }
    }
}
// 辅助 Kernel：将当前 Mesh 状态备份给中间变量，以便下一步使用
__global__ void copy_mesh_to_tmp(Element *d_Mesh, double d_U_tmp[Nx*Ny][NUM_VARS][N_MODE]) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < Nx * Ny) {
        for(int v=0; v<NUM_VARS; v++) 
            for(int k=0; k<N_MODE; k++)
                d_U_tmp[id][v][k] = d_Mesh[id].U[v][k];
    }
}

__global__ void apply_zhang_shu_limiter_kernel(Element *d_Mesh) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;

    Element *cell = &d_Mesh[id];
     double eps = 1e-13;
    double Ubar[NUM_VARS];
    for (int v = 0; v < NUM_VARS; v++) Ubar[v] = cell->U[v][0];

    if (Ubar[0] < eps) {
        for (int k = 1; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) cell->U[v][k] = 0.0;
        return;
    }

    double U_test[21][NUM_VARS];
    int pt = 0;
    for (int q = 0; q < N_QUAD; q++, pt++) {
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] += cell->U[v][k] * d_phi_vol[k][q];
    }
    for (int q = 0; q < 3; q++) {
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] += cell->U[v][k] * d_phi_face_T[k][q];
        pt++;
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] += cell->U[v][k] * d_phi_face_B[k][q];
        pt++;
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] += cell->U[v][k] * d_phi_face_L[k][q];
        pt++;
        for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] = 0.0;
        for (int k = 0; k < N_MODE; k++) for (int v = 0; v < NUM_VARS; v++) U_test[pt][v] += cell->U[v][k] * d_phi_face_R[k][q];
        pt++;
    }

    double rho_min = U_test[0][0];
    for (int i = 1; i < 21; i++) if (U_test[i][0] < rho_min) rho_min = U_test[i][0];

    double theta1 = 1.0;
    if (rho_min < eps) {
        theta1 = (Ubar[0] - eps) / (Ubar[0] - rho_min);
        if (theta1 < 0.0) theta1 = 0.0; if (theta1 > 1.0) theta1 = 1.0;
        for (int v = 0; v < NUM_VARS; v++) for (int k = 1; k < N_MODE; k++) cell->U[v][k] *= theta1;
        for (int i = 0; i < 21; i++) for (int v = 0; v < NUM_VARS; v++) U_test[i][v] = Ubar[v] + theta1 * (U_test[i][v] - Ubar[v]);
    }

    double theta2 = 1.0;
    for (int i = 0; i < 21; i++) {
        double p_i = calc_pressure_safe(U_test[i][0], U_test[i][1], U_test[i][2], U_test[i][3]);
        if (p_i < eps) {
            double t_L = 0.0, t_R = 1.0;
            for (int iter = 0; iter < 50; iter++) {
                double t_mid = 0.5 * (t_L + t_R);
                double Ut[NUM_VARS];
                for (int v = 0; v < NUM_VARS; v++) Ut[v] = Ubar[v] + t_mid * (U_test[i][v] - Ubar[v]);
                double p_mid = calc_pressure_safe(Ut[0], Ut[1], Ut[2], Ut[3]);
                if (p_mid < eps) t_R = t_mid; else t_L = t_mid;
            }
            if (t_L < theta2) theta2 = t_L;
        }
    }
    if (theta2 < 1.0)
        for (int v = 0; v < NUM_VARS; v++) for (int k = 1; k < N_MODE; k++) cell->U[v][k] *= theta2;
}

__global__ void copy_mesh_to_U0(Element *d_Mesh, double d_U0[Nx*Ny][NUM_VARS][N_MODE]) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < Nx * Ny) {
        for(int v=0; v<NUM_VARS; v++) for(int k=0; k<N_MODE; k++)
            d_U0[id][v][k] = d_Mesh[id].U[v][k];
    }
}

/* ============================================================
 * 辅助规约函子 (计算 dt)
 * ============================================================ */
struct WaveSpeedFunctor {
    Element* mesh;
    WaveSpeedFunctor(Element* _mesh) : mesh(_mesh) {}
    __device__ double operator()(const int& idx) const {
        double rho = mesh[idx].U[0][0], rhou = mesh[idx].U[1][0], rhov = mesh[idx].U[2][0], E = mesh[idx].U[3][0];
        double u = rhou / rho, v = rhov / rho;
        double p = calc_pressure(rho, rhou, rhov, E);
        double c = safe_sqrt(gamma_gas * p / rho);
   
        return fmax(fabs(u) + c, fabs(v) + c);
    }
};

struct DampSumFunctor {
    double (*damp)[3];
    DampSumFunctor(double (*_damp)[3]) : damp(_damp) {}
    __device__ double operator()(const int& idx) const {
        return damp[idx][0] + damp[idx][1] + damp[idx][2];
    }
};

/* ============================================================
 * CPU 端初始化与输出逻辑
 * ============================================================ */
static inline double legendre_1d(int i, double x) {
    if (i == 0) return 1.0; if (i == 1) return x; if (i == 2) return 0.5 * (3.0 * x * x - 1.0); return 0.0;
}
static inline double d_legendre_1d(int i, double x) {
    if (i == 0) return 0.0; if (i == 1) return 1.0; if (i == 2) return 3.0 * x; return 0.0;
}

void init_quadrature_and_matrices(void) {
    for (int j = 0; j < 3; j++) for (int i = 0; i < 3; i++) {
        int idx = j * 3 + i;
        r_quad_h[idx] = nodes1D_h[i]; s_quad_h[idx] = nodes1D_h[j]; w_quad_h[idx] = weights1D_h[i] * weights1D_h[j];
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
        for (int p = 0; p < 3; p++) {
            double np = nodes1D_h[p];
            phi_face_T_h[k][p] = legendre_1d(mk, np) * legendre_1d(nk, 1.0);
            phi_face_B_h[k][p] = legendre_1d(mk, np) * legendre_1d(nk, -1.0);
            phi_face_L_h[k][p] = legendre_1d(mk, -1.0) * legendre_1d(nk, np);
            phi_face_R_h[k][p] = legendre_1d(mk, 1.0) * legendre_1d(nk, np);
        }
    }
}
void init_condition(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh_h[jj * Nx + ii];
            
            for(int v = 0; v < NUM_VARS; v++) {
                for(int k = 0; k < N_MODE; k++) cell->U[v][k] = 0.0;
            }
            
            double xc = xL + (ii + 0.5) * dx;
            double yc = yL + (jj + 0.5) * dy;
            
            for (int q = 0; q < N_QUAD; q++) {
                double r = r_quad_h[q];
                double s = s_quad_h[q];
                double w = w_quad_h[q]; 
                
                double x_q = xc + r * (dx / 2.0);
                double y_q = yc + s * (dy / 2.0);
                
                // 【修改点】：植入目标平滑正弦波解析解 (t = 0)
                double rho = 1.0 + 0.2 * sin(M_PI * (x_q + y_q));
                double u = 0.7;
                double v = 0.3;
                double p = 1.0;
                
                double U_exact[NUM_VARS];
                U_exact[0] = rho;
                U_exact[1] = rho * u;
                U_exact[2] = rho * v;
                // 利用状态方程计算总能 E
                U_exact[3] = p / (gamma_gas - 1.0) + 0.5 * rho * (u * u + v * v);
                
                for (int k = 0; k < N_MODE; k++) {
                    double phi = phi_vol_h[k][q];
                    for (int var = 0; var < NUM_VARS; var++) {
                        cell->U[var][k] += w * U_exact[var] * phi;
                    }
                }
            }
            
            for (int k = 0; k < N_MODE; k++) {
                for (int var = 0; var < NUM_VARS; var++) {
                    cell->U[var][k] *= M_diag_inv_h[k];
                }
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
            int idx = jj * Nx + ii;
            Element *cell = &Mesh_h[idx]; 
            double xc = (ii + 0.5) * dx, yc = (jj + 0.5) * dy;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double r = nodes1D_h[j], s = nodes1D_h[i]; 
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
void calculate_error(double t) {
    double L1_err = 0.0;
    double L2_err = 0.0;
    double Linf_err = 0.0;
    double total_area = (xR - xL) * (yR - yL);

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh_h[jj * Nx + ii];
            double xc = xL + (ii + 0.5) * dx;
            double yc = yL + (jj + 0.5) * dy;

            for (int q = 0; q < N_QUAD; q++) {
                double r = r_quad_h[q];
                double s = s_quad_h[q];
                double w = w_quad_h[q]; 

                double x_q = xc + r * (dx / 2.0);
                double y_q = yc + s * (dy / 2.0);

                // 【修改点】：物理时间 t 下的精确解，波随时间向右上方移动
                double rho_exact = 1.0 + 0.2 * sin(M_PI * (x_q + y_q - t));

                double rho_num = 0.0;
                for (int k = 0; k < N_MODE; k++) {
                    rho_num += cell->U[0][k] * phi_vol_h[k][q];
                }

                double err = fabs(rho_num - rho_exact);
                double dV = (dx * dy / 4.0) * w;
                
                L1_err += err * dV;
                L2_err += (err * err) * dV;
                if (err > Linf_err) Linf_err = err;
            }
        }
    }

    L1_err = L1_err / total_area;
    L2_err = sqrt(L2_err / total_area);

    printf("\n================ Error Analysis ================\n");
    printf("Time       : %f\n", t);
    printf("L1 Error   : %e\n", L1_err);
    printf("L2 Error   : %e\n", L2_err);
    printf("Linf Error : %e\n", Linf_err);
    printf("================================================\n");
}
// ==== 新增: OEDG所需的全局标量提取函子 ====
struct VariableSumFunctor {
    Element* mesh;
    int var_idx;
    VariableSumFunctor(Element* _mesh, int _v) : mesh(_mesh), var_idx(_v) {}
    __device__ double operator()(const int& idx) const {
        return mesh[idx].U[var_idx][0]; 
    }
};

struct MaxDevFunctor {
    Element* mesh;
    double global_avg;
    int var_idx;
    MaxDevFunctor(Element* _mesh, double _avg, int _v) : mesh(_mesh), global_avg(_avg), var_idx(_v) {}
    __device__ double operator()(const int& idx) const {
        return fabs(mesh[idx].U[var_idx][0] - global_avg);
    }
};
int main(void) {
    init_quadrature_and_matrices();
    init_condition();

    // 拷贝常量到设备
    cudaMemcpyToSymbol(d_nodes1D, nodes1D_h, sizeof(nodes1D_h));
    cudaMemcpyToSymbol(d_weights1D, weights1D_h, sizeof(weights1D_h));
    cudaMemcpyToSymbol(d_mk_map, mk_map_h, sizeof(mk_map_h));
    cudaMemcpyToSymbol(d_nk_map, nk_map_h, sizeof(nk_map_h));
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

    // 分配 GPU 内存
    Element *d_Mesh, *d_Ghost_bottom, *d_Ghost_top, *d_Ghost_left, *d_Ghost_right;
    cudaMalloc(&d_Mesh, Nx * Ny * sizeof(Element));
    cudaMalloc(&d_Ghost_bottom, Nx * sizeof(Element)); cudaMalloc(&d_Ghost_top, Nx * sizeof(Element));
    cudaMalloc(&d_Ghost_left, Ny * sizeof(Element));   cudaMalloc(&d_Ghost_right, Ny * sizeof(Element));

    double (*d_cell_vertex_derivs)[4][6][NUM_VARS], (*d_damp_local)[3];
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
    
    int max_dim = (Nx > Ny) ? Nx : Ny;
    int gridMaxBoundary = (max_dim + blockSize1D - 1) / blockSize1D; 

    thrust::counting_iterator<int> iter(0);

    double current_time = 0.0; 
    int nit = 0;
    printf("%-10s  %-14s %-14s\n", "Step", "Time", "dt");
    printf("------------------------------------------\n");

while (current_time < T_END) {
        // 1. 仅计算边界和 Ghost cell 以获取正确的局部波速，用于计算准确的 dt
        compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
        apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);

        // 获取全场最大特征波速并计算时间步长
        double max_wave = thrust::transform_reduce(thrust::device, iter, iter + Nx * Ny, WaveSpeedFunctor(d_Mesh), 0.0, thrust::maximum<double>());
        if (max_wave < 1e-9) max_wave = 1.0;
        
        double dt = CFL * fmin(dx, dy) / max_wave;
        if (current_time + dt > T_END) dt = T_END - current_time;

        // 2. 将当前时间步的初始状态 U^n 备份到 d_U0
        copy_mesh_to_U0<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0);
        double max_dev_vars[NUM_VARS];

        // 3. 执行交替算子分裂的 SSP-RK3 推进
        for (int stage = 1; stage <= 3; stage++) {
            
            // (a) 利用 Thrust 提取全场的无量纲化分母 (全局平均值和最大偏差)
            for(int v = 0; v < NUM_VARS; v++) {
                double sum = thrust::transform_reduce(thrust::device, iter, iter + Nx*Ny, VariableSumFunctor(d_Mesh, v), 0.0, thrust::plus<double>());
                double global_avg = sum / (Nx * Ny);
                max_dev_vars[v] = thrust::transform_reduce(thrust::device, iter, iter + Nx*Ny, MaxDevFunctor(d_Mesh, global_avg, v), 0.0, thrust::maximum<double>());
            }

            // (b) 备份中间状态: 如果是 Stage 2 或 3，当前的 d_Mesh 是上一阶段的输出，需备份至 d_U_tmp 参与组合推进
            if (stage > 1) {
                copy_mesh_to_tmp<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U_tmp); 
            }

            // (c) 更新当前状态的边界信息与单元顶点导数
            compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
            apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
            precompute_vertex_derivs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_cell_vertex_derivs);
            
            // (d) 计算基于全场尺度的 OEDG 阻尼系数
            compute_oedg_damp_kernel<<<gridSize2D, blockSize2D>>>(
                d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, 
                d_cell_vertex_derivs, d_damp_local, 
                max_dev_vars[0], max_dev_vars[1], max_dev_vars[2], max_dev_vars[3]
            );

            // (e) 评估双曲空间通量右端项 L_f(U)
            compute_rhs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_RHS);

            // (f) 纯物理推进: 根据当前的 stage，严格按照 SSP-RK3 组装状态 U^*
            // 当 stage == 1 时，d_U_tmp 未使用，传入 d_U0 占位即可
            rk3_physical_stage_kernel<<<gridSize1D, blockSize1D>>>(
                d_Mesh, d_U0, (stage == 1 ? d_U0 : d_U_tmp), d_RHS, stage, dt
            );

            // (g) 纯数学滤波: 将 OE 精确算子施加于物理推进后的中间态 U^*，得到无振荡状态，0阶模态天然守恒
            apply_oedg_filter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_damp_local, dt);
            
            // (h) 施加保正限制器以维持极低密度的物理有效性
            apply_zhang_shu_limiter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh);
        }

        current_time += dt;
        nit++;
        printf("%-10d  %-14.6e %-14.6e \n", nit, current_time, dt);
    }
    cudaMemcpy(Mesh_h, d_Mesh, Nx * Ny * sizeof(Element), cudaMemcpyDeviceToHost);
    
    printf("GPU Computation Complete.\n");
    output_results(T_END);
  calculate_error(T_END);
    cudaFree(d_Mesh); cudaFree(d_Ghost_bottom); cudaFree(d_Ghost_top); cudaFree(d_Ghost_left); cudaFree(d_Ghost_right);
    cudaFree(d_cell_vertex_derivs); cudaFree(d_damp_local);
    cudaFree(d_U0); cudaFree(d_U_tmp);
    cudaFree(d_RHS);
    return 0;
}