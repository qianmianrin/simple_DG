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

#define Nx 320
#define Ny 160
#define xL 0.0
#define xR 1.0
#define yL -0.25
#define yR 0.25
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma_gas (5.0 / 3.0)
#define NUM_VARS 4 
#define Mach 800.0 


#define CFL   0.01
#define T_END 0.001


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
    if (x < -1e-8) return sqrt(x); 
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

__device__ static void build_eigen_x(double rho, double u, double v, double p, double L[NUM_VARS][NUM_VARS]) {
    double c = sqrt(gamma_gas * p / rho);
    double B1 = gamma_gas - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    double f_L = B1 / c;

    // 严格匹配论文 Eq (2.14), 法向为 X (n1=1, n2=0)
    L[0][0] = f_L * 0.5 * (B2 + u * c);
    L[0][1] = f_L * -0.5 * (B1 * u + c);
    L[0][2] = f_L * -0.5 * (B1 * v);
    L[0][3] = f_L * 0.5 * B1;

    L[1][0] = f_L * (c * c - B2);
    L[1][1] = f_L * B1 * u;
    L[1][2] = f_L * B1 * v;
    L[1][3] = f_L * -B1;

    L[2][0] = f_L * (-v * c);
    L[2][1] = 0.0;
    L[2][2] = f_L * c;
    L[2][3] = 0.0;

    L[3][0] = f_L * 0.5 * (B2 - u * c);
    L[3][1] = f_L * -0.5 * (B1 * u - c);
    L[3][2] = f_L * -0.5 * (B1 * v);
    L[3][3] = f_L * 0.5 * B1;
}

__device__ static void build_eigen_y(double rho, double u, double v, double p, double L[NUM_VARS][NUM_VARS]) {
    double c = sqrt(gamma_gas * p / rho);
    double B1 = gamma_gas - 1.0;
    double B2 = 0.5 * B1 * (u * u + v * v);
    double f_L = B1 / c;

    // 严格匹配论文 Eq (2.14), 法向为 Y (n1=0, n2=1)
    L[0][0] = f_L * 0.5 * (B2 + v * c);
    L[0][1] = f_L * -0.5 * (B1 * u);
    L[0][2] = f_L * -0.5 * (B1 * v + c);
    L[0][3] = f_L * 0.5 * B1;

    L[1][0] = f_L * (c * c - B2);
    L[1][1] = f_L * B1 * u;
    L[1][2] = f_L * B1 * v;
    L[1][3] = f_L * -B1;

    L[2][0] = f_L * (u * c);
    L[2][1] = f_L * (-c);
    L[2][2] = 0.0;
    L[2][3] = 0.0;

    L[3][0] = f_L * 0.5 * (B2 - v * c);
    L[3][1] = f_L * -0.5 * (B1 * u);
    L[3][2] = f_L * -0.5 * (B1 * v - c);
    L[3][3] = f_L * 0.5 * B1;
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

__device__ static void get_primitive_vars(double U[NUM_VARS], double *rho, double *u, double *v, double *p) {
    *rho = U[0]; *u = U[1] / U[0]; *v = U[2] / U[0];
    *p = calc_pressure(U[0], U[1], U[2], U[3]);
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
    
    double rho1 = 5.0, u1 = Mach, v1 = 0.0, p1 = 0.4127;
    double E1 = p1 / (gamma_gas - 1.0) + 0.5 * rho1 * (u1 * u1 + v1 * v1);
    double rho0 = 0.5, u0 = 0.0, v0 = 0.0, p0 = 0.4127;
    double E0 = p0 / (gamma_gas - 1.0);

    if (id < Nx) {
        for(int v=0; v<NUM_VARS; v++) for(int k=0; k<N_MODE; k++) {
            d_Ghost_bottom[id].U[v][k] = d_Mesh[0 * Nx + id].U[v][k];
            d_Ghost_top[id].U[v][k] = d_Mesh[(Ny - 1) * Nx + id].U[v][k];
        }
        compute_face_values_device(&d_Ghost_bottom[id]);
        compute_face_values_device(&d_Ghost_top[id]);
    }
    
    if (id < Ny) {
        for(int v=0; v<NUM_VARS; v++) for(int k=0; k<N_MODE; k++) {
            d_Ghost_right[id].U[v][k] = d_Mesh[id * Nx + (Nx - 1)].U[v][k];
            d_Ghost_left[id].U[v][k] = 0.0;
        }
        double y_center = yL + (id + 0.5) * dy;
        if (y_center >= -0.05 && y_center <= 0.05) {
            d_Ghost_left[id].U[0][0] = rho1; d_Ghost_left[id].U[1][0] = rho1 * u1;
            d_Ghost_left[id].U[2][0] = rho1 * v1; d_Ghost_left[id].U[3][0] = E1;
        } else {
            d_Ghost_left[id].U[0][0] = rho0; d_Ghost_left[id].U[1][0] = rho0 * u0;
            d_Ghost_left[id].U[2][0] = rho0 * v0; d_Ghost_left[id].U[3][0] = E0;
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

__global__ void compute_damp_coeffs_kernel(Element *d_Mesh, Element *d_Ghost_bottom, Element *d_Ghost_top, 
                                           Element *d_Ghost_left, Element *d_Ghost_right,
                                           double d_cell_vertex_derivs[Nx*Ny][4][6][NUM_VARS],
                                           double d_damp_local[Nx*Ny][3]) {
    int ii = blockIdx.x * blockDim.x + threadIdx.x;
    int jj = blockIdx.y * blockDim.y + threadIdx.y;
    if (ii >= Nx || jj >= Ny) return;
    int idx = jj * Nx + ii;

    Element *C_cell = &d_Mesh[idx];
    Element *B_cell = (jj == 0) ? &d_Ghost_bottom[ii] : &d_Mesh[(jj - 1) * Nx + ii];
    Element *T_cell = (jj == Ny - 1) ? &d_Ghost_top[ii] : &d_Mesh[(jj + 1) * Nx + ii];
    Element *L_cell = (ii == 0) ? &d_Ghost_left[jj] : &d_Mesh[jj * Nx + (ii - 1)];
    Element *R_cell = (ii == Nx - 1) ? &d_Ghost_right[jj] : &d_Mesh[jj * Nx + (ii + 1)];

    double deriv_curr[4][6][NUM_VARS]; 
    double deriv_B[2][6][NUM_VARS], deriv_T[2][6][NUM_VARS];    
    double deriv_L[2][6][NUM_VARS], deriv_R[2][6][NUM_VARS];    

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

    double L_mats[4][2][NUM_VARS][NUM_VARS];
    for (int vert = 0; vert < 4; vert++) {
        for (int edge = 0; edge < 2; edge++) {
            double U_mean_curr[NUM_VARS], U_mean_neb[NUM_VARS];
            for (int v = 0; v < NUM_VARS; v++) {
                U_mean_curr[v] = C_cell->U[v][0];
                if (edge == 0) { 
                    if (vert == 0 || vert == 3) U_mean_neb[v] = L_cell->U[v][0]; else U_mean_neb[v] = R_cell->U[v][0];
                } else { 
                    if (vert == 0 || vert == 1) U_mean_neb[v] = B_cell->U[v][0]; else U_mean_neb[v] = T_cell->U[v][0];
                }
            }
            double rho_c, u_c, v_c, p_c, rho_n, u_n, v_n, p_n;
            get_primitive_vars(U_mean_curr, &rho_c, &u_c, &v_c, &p_c);
            get_primitive_vars(U_mean_neb, &rho_n, &u_n, &v_n, &p_n);

            double H_c = (U_mean_curr[3] + p_c) / rho_c; double H_n = (U_mean_neb[3] + p_n) / rho_n;
            double R = safe_sqrt(rho_n / rho_c); double inv_R_1 = 1.0 / (R + 1.0);
            
            double rho_roe = R * rho_c; 
            double u_roe = (u_c + R * u_n) * inv_R_1; double v_roe = (v_c + R * v_n) * inv_R_1;
            double H_roe = (H_c + R * H_n) * inv_R_1;
            
            double vel_sq = u_roe * u_roe + v_roe * v_roe;
            double p_roe = (gamma_gas - 1.0) * (H_roe - 0.5 * vel_sq) * rho_roe / gamma_gas;
            
            if (edge == 0) build_eigen_x(rho_roe, u_roe, v_roe, p_roe, L_mats[vert][edge]);
            else           build_eigen_y(rho_roe, u_roe, v_roe, p_roe, L_mats[vert][edge]);
        }
    }

    int order_map[6] = {0, 1, 1, 2, 2, 2}; 
    for (int l = 0; l <= 2; l++) {
        double max_wave_jump = 0.0;
        for (int s = 0; s < NUM_VARS; s++) {
            double sum_alpha = 0.0; 
            for (int d = 0; d < 6; d++) {
                if (order_map[d] != l) continue;
                double sum_vertices = 0.0; 
                for (int vert = 0; vert < 4; vert++) {
                    double v_jump_sq = 0.0;
                    for (int edge = 0; edge < 2; edge++) {
                        double c_jump = 0.0;
                        for (int v_idx = 0; v_idx < NUM_VARS; v_idx++) c_jump += L_mats[vert][edge][s][v_idx] * jump[vert][edge][d][v_idx];
                        v_jump_sq += c_jump * c_jump;
                    }
                    sum_vertices += v_jump_sq;
                }
                sum_alpha += safe_sqrt(sum_vertices / 4.0);
            }
            if (sum_alpha > max_wave_jump) max_wave_jump = sum_alpha;
        }
        
        double factorial = (l==0)? 1.0 : ((l==1)? 1.0 : 2.0);
        double h_pow = (l==0)? 1.0 : ((l==1)? dx : dx*dx);
        d_damp_local[idx][l] = ((2.0 * (2 * l + 1)) / 3.0) * h_pow / factorial * max_wave_jump; 
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

// 依据论文 Eq (2.24) 实现的指数型 RK3 积分器
// 彻底取代了原先将阻尼项作为显式 RHS 的不当做法
__global__ void update_exp_rk3_stage_kernel(Element *d_Mesh, 
                                            double d_U0[Nx*Ny][NUM_VARS][N_MODE], 
                                            double d_U_prev[Nx*Ny][NUM_VARS][N_MODE], 
                                            double d_RHS[Nx*Ny][NUM_VARS][N_MODE], 
                                            double d_damp_local[Nx*Ny][3],
                                            int stage, double dt) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= Nx * Ny) return;
    
    double h_K = fmax(dx, dy);
    double damp0 = d_damp_local[id][0];
    double damp1 = d_damp_local[id][1];
    double damp2 = d_damp_local[id][2];
    
    // 论文中取 a0 为当前单元内各阶阻尼的累加最大值尺度
    double a0 = (damp0 + damp1 + damp2) / h_K;
    double z = a0 * dt;
    
    // 泰勒展开近似指数积分因子，避免直接计算 exp() 带来的浮点灾难
    double s1 = 1.0 + z + 0.5 * z * z + (1.0 / 6.0) * z * z * z;
    double s2 = 1.0 + 0.5 * z + 0.125 * z * z + (1.0 / 48.0) * z * z * z;

    int order_map[6] = {0, 1, 1, 2, 2, 2};
    double sigma_k[N_MODE];
    sigma_k[0] = 0.0; // 守恒律核心：单元平均值(0阶模态)绝对不加任何阻尼
    for(int k = 1; k < N_MODE; k++) {
        double sum_sigma = 0.0;
        for(int l = 0; l <= order_map[k]; l++) {
            if(l == 0) sum_sigma += damp0 / h_K;
            if(l == 1) sum_sigma += damp1 / h_K;
            if(l == 2) sum_sigma += damp2 / h_K;
        }
        sigma_k[k] = sum_sigma;
    }

    for (int v = 0; v < NUM_VARS; v++) {
        // 0阶模态（守恒量平均值）走标准 RK3，不受阻尼影响
        if (stage == 1) {
            d_Mesh[id].U[v][0] = d_U0[id][v][0] + dt * d_RHS[id][v][0];
            for (int k = 1; k < N_MODE; k++) {
                double rhs_mod = d_RHS[id][v][k] - sigma_k[k] * d_U0[id][v][k];
                d_Mesh[id].U[v][k] = (1.0 / s1) * (d_U0[id][v][k] + dt * (rhs_mod + a0 * d_U0[id][v][k]));
            }
        } else if (stage == 2) {
            d_Mesh[id].U[v][0] = 0.75 * d_U0[id][v][0] + 0.25 * (d_U_prev[id][v][0] + dt * d_RHS[id][v][0]);
            for (int k = 1; k < N_MODE; k++) {
                double rhs_mod = d_RHS[id][v][k] - sigma_k[k] * d_U_prev[id][v][k];
                d_Mesh[id].U[v][k] = (0.75 / s2) * d_U0[id][v][k] + (s1 / (4.0 * s2)) * (d_U_prev[id][v][k] + dt * (rhs_mod + a0 * d_U_prev[id][v][k]));
            }
        } else if (stage == 3) {
            d_Mesh[id].U[v][0] = (1.0 / 3.0) * d_U0[id][v][0] + (2.0 / 3.0) * (d_U_prev[id][v][0] + dt * d_RHS[id][v][0]);
            for (int k = 1; k < N_MODE; k++) {
                double rhs_mod = d_RHS[id][v][k] - sigma_k[k] * d_U_prev[id][v][k];
                d_Mesh[id].U[v][k] = (1.0 / (3.0 * s1)) * d_U0[id][v][k] + (2.0 * s2 / (3.0 * s1)) * (d_U_prev[id][v][k] + dt * (rhs_mod + a0 * d_U_prev[id][v][k]));
            }
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
    const double eps = 1e-13;
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
        double c_jet = safe_sqrt(gamma_gas * 0.4127 / 5.0);
        return fmax(fmax(fabs(u) + c, Mach + c_jet), fabs(v) + c);
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
            for(int v = 0; v < NUM_VARS; v++) for(int k = 0; k < N_MODE; k++) cell->U[v][k] = 0.0;
            double rho = 0.5, u = 0.0, v = 0.0, p = 0.4127;
            cell->U[0][0] = rho; cell->U[1][0] = rho * u; cell->U[2][0] = rho * v;
            cell->U[3][0] = p / (gamma_gas - 1.0) + 0.5 * rho * (u * u + v * v);
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
        // 先计算边界和 Ghost cell 以计算准确的 dt
        compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
        apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
        precompute_vertex_derivs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_cell_vertex_derivs);
        compute_damp_coeffs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs, d_damp_local);

        // 获取全场最大特征波速
        double max_wave = thrust::transform_reduce(thrust::device, iter, iter + Nx * Ny, WaveSpeedFunctor(d_Mesh), 0.0, thrust::maximum<double>());
        if (max_wave < 1e-9) max_wave = 1.0;
        
        // 移除了 max_a0 对 dt 的限制约束，仅保留声速物理约束
        double dt = CFL * fmin(dx, dy) / max_wave;
        if (current_time + dt > T_END) dt = T_END - current_time;

        copy_mesh_to_U0<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0);

        // ==========================================
        // Stage 1: U^{(1)} 
        // 评估时间: current_time
        // ==========================================
        compute_rhs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_RHS);
        update_exp_rk3_stage_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0, d_U0, d_RHS, d_damp_local, 1, dt);
        apply_zhang_shu_limiter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh);

        // ==========================================
        // Stage 2: U^{(2)}
        // 评估时间: current_time + dt
        // ==========================================
        copy_mesh_to_tmp<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U_tmp); 
        
        compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
        apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
        precompute_vertex_derivs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_cell_vertex_derivs);
        compute_damp_coeffs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs, d_damp_local);
        
        compute_rhs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_RHS);
        update_exp_rk3_stage_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0, d_U_tmp, d_RHS, d_damp_local, 2, dt);
        apply_zhang_shu_limiter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh);

        // ==========================================
        // Stage 3: U^{(n+1)}
        // 评估时间: current_time + 0.5 * dt
        // ==========================================
        copy_mesh_to_tmp<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U_tmp);
        
        compute_boundary_faces_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh);
        apply_ghost_cells_kernel<<<gridMaxBoundary, blockSize1D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right);
        precompute_vertex_derivs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_cell_vertex_derivs);
        compute_damp_coeffs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_cell_vertex_derivs, d_damp_local);
        
        compute_rhs_kernel<<<gridSize2D, blockSize2D>>>(d_Mesh, d_Ghost_bottom, d_Ghost_top, d_Ghost_left, d_Ghost_right, d_RHS);
        update_exp_rk3_stage_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh, d_U0, d_U_tmp, d_RHS, d_damp_local, 3, dt);
        apply_zhang_shu_limiter_kernel<<<gridSize1D, blockSize1D>>>(d_Mesh);

        current_time += dt; 
        nit++;
        printf("%-10d  %-14.6e %-14.6e \n", nit, current_time, dt);
    }

    cudaMemcpy(Mesh_h, d_Mesh, Nx * Ny * sizeof(Element), cudaMemcpyDeviceToHost);
    
    printf("GPU Computation Complete.\n");
    output_results(T_END);

    cudaFree(d_Mesh); cudaFree(d_Ghost_bottom); cudaFree(d_Ghost_top); cudaFree(d_Ghost_left); cudaFree(d_Ghost_right);
    cudaFree(d_cell_vertex_derivs); cudaFree(d_damp_local);
    cudaFree(d_U0); cudaFree(d_U_tmp);
    cudaFree(d_RHS);
    return 0;
}