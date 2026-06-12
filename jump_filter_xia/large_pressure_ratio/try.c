#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define GAMMA 1.4
#define NVAR 3
#define lambdacoff 1
 
typedef struct {
    int Nx;
    double *U_ext;     // (Nx+6) * NVAR
    double *Fhat;      // (Nx+1) * NVAR
    
    // 用于 1D WENO 计算的临时缓冲区
    double *ws_double; 
} SolverWorkspace;

// 申请工作空间
SolverWorkspace* create_workspace(int Nx) {
    SolverWorkspace *ws = (SolverWorkspace*)malloc(sizeof(SolverWorkspace));
    ws->Nx = Nx;
    ws->U_ext = (double*)malloc((Nx + 6) * NVAR * sizeof(double));
    ws->Fhat  = (double*)malloc((Nx + 1) * NVAR * sizeof(double));
    
    int max_dim = Nx + 6;
    // 分配足够的一维临时数组: rho_avg, u_avg, H_avg, c_avg (各 Nf), lambda (3*Nf) 
    ws->ws_double = (double*)malloc(max_dim * 10 * sizeof(double)); 
    return ws;
}

void free_workspace(SolverWorkspace *ws) {
    free(ws->U_ext);
    free(ws->Fhat);
    free(ws->ws_double);
    free(ws);
}

static inline double pressure(const double *w) {
    double rho = w[0];
    double mx  = w[1];
    double E   = w[2];
    return (GAMMA - 1.0) * (E - 0.5 * (mx * mx) / rho);
}

static inline double solve_t(const double *h, const double *w_j, double eps) {
    const double gm1 = GAMMA - 1.0;          
    double rho0  = h[0], rhou0 = h[1], E0 = h[2];
    double drho  = w_j[0] - rho0, drhou = w_j[1] - rhou0, dE = w_j[2] - E0;
    double coef  = 2.0 * eps / gm1;            

    double a =  2.0*drho*dE - drhou*drhou;
    double b =  2.0*(rho0*dE + E0*drho - rhou0*drhou) - coef * drho;
    double c =  2.0*rho0*E0 - rhou0*rhou0 - coef * rho0;

    if (c < 0.0) return 0.0; 

    if (fabs(a) < 1e-14) {
        if (fabs(b) < 1e-14) return 1.0;
        double t = -c / b;
        return (t >= 0.0 && t <= 1.0) ? t : ((t < 0.0) ? 0.0 : 1.0);
    }

    double disc = b*b - 4.0*a*c;
    if (disc < 0.0) return 1.0;

    double sq = sqrt(disc);
    double q_root = -0.5 * (b + (b > 0.0 ? sq : -sq));
    double t1 = q_root / a;
    double t2 = c / q_root;

    if (t1 > t2) { double tmp=t1; t1=t2; t2=tmp; }

    double t_sol = (a > 0.0) ? t1 : t2;

    if (t_sol < 0.0) t_sol = 0.0;
    if (t_sol > 1.0) t_sol = 1.0;

    return t_sol;
}

static inline void pp_limiter(double *face, double *h, double *target, double eps) {
    double omega = 1.0/12.0;
    double q[3] = {0.0};
    for(int m=0; m<3; m++){
        q[m] = (1.0/(1.0-omega)) * (h[m] - omega*face[m]);
    }
    double rho_min = fmin(q[0], face[0]);

    double theta1 = 1.0;
    if (rho_min < eps) {
        theta1 = (h[0] - eps) / (h[0] - rho_min); 
    }

    double rho_star = theta1 * (face[0] - h[0]) + h[0];
    double face_star[3]; 
    face_star[0] = rho_star;
    face_star[1] = face[1];
    face_star[2] = face[2];

    double q_star[3] = {0.0};
    for (int m = 0; m < 3; m++) {
        q_star[m] = (1.0/(1.0-omega)) * (h[m] - omega*face_star[m]);
    }

    double t1 = 1.0, t2 = 1.0;
    if(pressure(face_star) < 1e-13){
        t1 = solve_t(h, face_star, 1e-13);
    }

    if(pressure(q_star) < 1e-13){
        t2 = solve_t(h, q_star, 1e-13);
    }

    double theta2 = fmin(t1, t2);


    for(int m=0; m<3; m++){
        target[m] = theta2 * (face_star[m] - h[m]) + h[m];
    }
}

void apply_BC_1D(double *U, SolverWorkspace *ws) {
    int Nx = ws->Nx;
    double *U_ext = ws->U_ext; 

    // 内部场拷贝
    for (int i = 0; i < Nx; i++) {
        memcpy(&U_ext[(i + 3) * NVAR], &U[i * NVAR], NVAR * sizeof(double));
    }

    // 零阶外推边界 (Transmissive)
    for (int k = 0; k < NVAR; k++) {
        // 左边界
        U_ext[0 * NVAR + k] = U[0 * NVAR + k];
        U_ext[1 * NVAR + k] = U[0 * NVAR + k];
        U_ext[2 * NVAR + k] = U[0 * NVAR + k];
        
        // 右边界
        U_ext[(Nx + 3) * NVAR + k] = U[(Nx - 1) * NVAR + k];
        U_ext[(Nx + 4) * NVAR + k] = U[(Nx - 1) * NVAR + k];
        U_ext[(Nx + 5) * NVAR + k] = U[(Nx - 1) * NVAR + k];
    }
}

static inline double weno5_z(double v0, double v1, double v2, double v3, double v4, int is_right)
{
    const double eps = 1e-12;
    double f0, f1, f2;
    double vv0,vv1,vv2,vv3,vv4;
    vv0=v0,vv1=v1,vv2=v2,vv3=v3,vv4=v4;
    double vmax=v0;
    double vmin=v0;
  
if (vmax<v0) vmax=v0;
if (vmax<v1) vmax=v1;
   if (vmax<v2) vmax=v2;
if (vmax<v3) vmax=v3;
if (vmax<v4) vmax=v4;
if (vmin>v0) vmin=v0;
if (vmin > v1) vmin = v1;
if (vmin > v2) vmin = v2;
if (vmin > v3) vmin = v3;
if (vmin > v4) vmin = v4;
double diff = vmax - vmin;
if (diff < 1e-9) { // 设一个极小的阈值 epsilon
    vv0 = vv1 = vv2 = vv3 = vv4 = 0.0; // 或者根据业务逻辑设
} else {
vv0=(v0-vmin)/(vmax-vmin);
vv1 = (v1 - vmin) / (vmax - vmin);
vv2 = (v2 - vmin) / (vmax - vmin);
vv3 = (v3 - vmin) / (vmax - vmin);
vv4 = (v4 - vmin) / (vmax - vmin);
}
    // 1. Reconstruction polynomials (保持原样)
    if (is_right) {
        f0 = (2.0*v0 - 7.0*v1 + 11.0*v2) * 0.1666666666666667;
        f1 = (-v1 + 5.0*v2 + 2.0*v3) * 0.1666666666666667;
        f2 = (2.0*v2 + 5.0*v3 - v4) * 0.1666666666666667;
    } else {
        f0 = (-v0 + 5.0*v1 + 2.0*v2) * 0.1666666666666667;
        f1 = (2.0*v1 + 5.0*v2 - v3) * 0.1666666666666667;
        f2 = (11.0*v2 - 7.0*v3 + 2.0*v4) * 0.1666666666666667;
    }

    // 2. Smoothness Indicators (IS) (保持原样)
    double d1 = vv0 - 2.0*vv1 + vv2;
    double d2 = vv0 - 4.0*vv1 + 3.0*vv2;
    double beta0 = 1.0833333333333333 * d1*d1 + 0.25 * d2*d2;

    d1 = vv1 - 2.0*vv2 + vv3;
    d2 = vv1 - vv3;
    double beta1 = 1.0833333333333333 * d1*d1 + 0.25 * d2*d2;

    d1 = vv2 - 2.0*vv3 + vv4;
    d2 = 3.0*vv2 - 4.0*vv3 + vv4;
    double beta2 = 1.0833333333333333 * d1*d1 + 0.25 * d2*d2;

    // 3. WENO-Z Nonlinear weights
    double ideal_w0 = is_right ? 0.1 : 0.3;
    double ideal_w1 = 0.6;
    double ideal_w2 = is_right ? 0.3 : 0.1;

    // WENO-Z 核心：计算高阶光滑度指示子 tau5 = |beta0 - beta2|
    double tau5 = fabs(beta0 - beta2);

    double s0 = eps + beta0;
    double s1 = eps + beta1;
    double s2 = eps + beta2;

    // 计算 WENO-Z 权重: alpha = d * (1 + (tau5/s)^p), 此处采用标准 p=2
    // 相比 JS-WENO 的 d/(s^2)，WENO-Z 的分母永远不会趋近于 0，计算更稳定
    double alpha0 = ideal_w0 * (1.0 + (tau5 * tau5) / (s0 * s0));
    double alpha1 = ideal_w1 * (1.0 + (tau5 * tau5) / (s1 * s1));
    double alpha2 = ideal_w2 * (1.0 + (tau5 * tau5) / (s2 * s2));
    
    // 4. Combination
    return (alpha0*f0 + alpha1*f1 + alpha2*f2) / (alpha0 + alpha1 + alpha2);
}



void compute_Fhat_1D(SolverWorkspace *ws) {
    int Nx = ws->Nx;
    int Nf = Nx + 1;
    double *U_ext = ws->U_ext;
    
    double *ptr = ws->ws_double;
    double *u_avg = ptr; ptr += Nf;
    double *H_avg = ptr; ptr += Nf;
    double *c_avg = ptr; ptr += Nf;

    for(int i = 0; i < Nf; i++) {
        int idxL = i + 2;
        int idxR = i + 3;

        double *UL = &U_ext[idxL * NVAR];
        double *UR = &U_ext[idxR * NVAR];

        double rhoL = UL[0], rhoR = UR[0];
        double inv_rhoL = 1.0/rhoL, inv_rhoR = 1.0/rhoR;
        
        double uL = UL[1] * inv_rhoL;
        double uR = UR[1] * inv_rhoR;

        double pL = (GAMMA - 1.0) * (UL[2] - 0.5 * rhoL * uL*uL);
        double pR = (GAMMA - 1.0) * (UR[2] - 0.5 * rhoR * uR*uR);
    
        double HL = (UL[2] + pL) * inv_rhoL;
        double HR = (UR[2] + pR) * inv_rhoR;

        double srL = sqrt(rhoL);
        double srR = sqrt(rhoR);
        double inv_sr_sum = 1.0 / (srL + srR);

        u_avg[i] = (srL * uL + srR * uR) * inv_sr_sum;
        H_avg[i] = (srL * HL + srR * HR) * inv_sr_sum;
        c_avg[i] = sqrt((GAMMA - 1.0) * (H_avg[i] - 0.5 * u_avg[i]*u_avg[i]));
    }

    double alpha_global = 1e-8;
    for (int i = 0; i < Nx + 6; i++) {
        double *U = &U_ext[i * NVAR];
        double rho = U[0];
        double u = U[1] / rho;
        double p = (GAMMA - 1.0) * (U[2] - 0.5 * rho * u*u);
        double c = sqrt(GAMMA * p / rho);
        double lam = fabs(u) + c; 
        if (lam > alpha_global) alpha_global = lam;
    }

    double eps1 = 1e-13, eps2 = 1e-13;
    for(int i = 0; i < Nf; i++) {
        eps1 = fmin(eps1, U_ext[(i+2)*NVAR + 0]); 
        eps2 = fmin(eps2, U_ext[(i+3)*NVAR + 0]);
    }

    for(int i = 0; i < Nf; i++) {
        int iL = i + 2; 
        double alpha = alpha_global;
        double u = u_avg[i], H = H_avg[i], c = c_avg[i];
        
        double rc = 1.0 / c;
        double b1 = (GAMMA - 1.0) * rc * rc;
        double b2 = 0.5 * u*u * b1;
        double t1 = b1 * u;
        
        // 3x3 L matrix
        double L00 = 0.5 * (b2 + u*rc), L01 = -0.5 * (t1 + rc), L02 = 0.5 * b1;
        double L10 = 1.0 - b2,          L11 = t1,               L12 = -b1;
        double L20 = 0.5 * (b2 - u*rc), L21 = -0.5 * (t1 - rc), L22 = 0.5 * b1;

        double F_interface_plus[3] = {0.0};
        double F_interface_minus[3] = {0.0};

        for(int m = 0; m < 3; m++) {
            double fp_stencil[5], fm_stencil[5];

            for(int s = 0; s < 5; s++) {
                int idx = iL - 2 + s;
                double *U = &U_ext[idx*NVAR];
                double phys_F[3] = { U[1], U[1]*U[1]/U[0] + (GAMMA-1.0)*(U[2]-0.5*U[1]*U[1]/U[0]), U[1]/U[0]*(U[2] + (GAMMA-1.0)*(U[2]-0.5*U[1]*U[1]/U[0])) };
                
                double w_val = 0.0, fc_val = 0.0;
                if (m == 0) { w_val = L00*U[0] + L01*U[1] + L02*U[2]; fc_val = L00*phys_F[0] + L01*phys_F[1] + L02*phys_F[2]; } 
                else if (m == 1) { w_val = L10*U[0] + L11*U[1] + L12*U[2]; fc_val = L10*phys_F[0] + L11*phys_F[1] + L12*phys_F[2]; } 
                else { w_val = L20*U[0] + L21*U[1] + L22*U[2]; fc_val = L20*phys_F[0] + L21*phys_F[1] + L22*phys_F[2]; }
                fp_stencil[s] = 0.5 * (w_val + fc_val/alpha);
            }
            
            for(int s = 0; s < 5; s++) {
                int idx = iL - 1 + s;
                double *U = &U_ext[idx*NVAR];
                double phys_F[3] = { U[1], U[1]*U[1]/U[0] + (GAMMA-1.0)*(U[2]-0.5*U[1]*U[1]/U[0]), U[1]/U[0]*(U[2] + (GAMMA-1.0)*(U[2]-0.5*U[1]*U[1]/U[0])) };
                 
                double w_val = 0.0, fc_val = 0.0;
                if (m == 0) { w_val = L00*U[0] + L01*U[1] + L02*U[2]; fc_val = L00*phys_F[0] + L01*phys_F[1] + L02*phys_F[2]; } 
                else if (m == 1) { w_val = L10*U[0] + L11*U[1] + L12*U[2]; fc_val = L10*phys_F[0] + L11*phys_F[1] + L12*phys_F[2]; } 
                else { w_val = L20*U[0] + L21*U[1] + L22*U[2]; fc_val = L20*phys_F[0] + L21*phys_F[1] + L22*phys_F[2]; }
                fm_stencil[s] = 0.5 * (w_val - fc_val/alpha);
            }

            double f_plus = weno5_z(fp_stencil[0], fp_stencil[1], fp_stencil[2], fp_stencil[3], fp_stencil[4], 1);
            double f_minus = weno5_z(fm_stencil[0], fm_stencil[1], fm_stencil[2], fm_stencil[3], fm_stencil[4], 0);

            // Project back using R matrix
            if (m == 0) { 
                F_interface_plus[0] += f_plus;          F_interface_minus[0] += f_minus;
                F_interface_plus[1] += f_plus*(u-c);    F_interface_minus[1] += f_minus*(u-c);
                F_interface_plus[2] += f_plus*(H-u*c);  F_interface_minus[2] += f_minus*(H-u*c);
            } else if (m == 1) { 
                F_interface_plus[0] += f_plus;          F_interface_minus[0] += f_minus;
                F_interface_plus[1] += f_plus*u;        F_interface_minus[1] += f_minus*u;
                F_interface_plus[2] += f_plus*0.5*u*u;  F_interface_minus[2] += f_minus*0.5*u*u;
            } else { 
                F_interface_plus[0] += f_plus;          F_interface_minus[0] += f_minus;
                F_interface_plus[1] += f_plus*(u+c);    F_interface_minus[1] += f_minus*(u+c);
                F_interface_plus[2] += f_plus*(H+u*c);  F_interface_minus[2] += f_minus*(H+u*c);
            }
        } 

        double F_interface_plus_re[3] = {0.0}, F_interface_minus_re[3] = {0.0};
        double *UL = &U_ext[iL*NVAR];
        double phys_FL[3] = { UL[1], UL[1]*UL[1]/UL[0] + (GAMMA-1.0)*(UL[2]-0.5*UL[1]*UL[1]/UL[0]), UL[1]/UL[0]*(UL[2] + (GAMMA-1.0)*(UL[2]-0.5*UL[1]*UL[1]/UL[0])) };
        
        double *UR = &U_ext[(iL+1)*NVAR];
        double phys_FR[3] = { UR[1], UR[1]*UR[1]/UR[0] + (GAMMA-1.0)*(UR[2]-0.5*UR[1]*UR[1]/UR[0]), UR[1]/UR[0]*(UR[2] + (GAMMA-1.0)*(UR[2]-0.5*UR[1]*UR[1]/UR[0])) };

        double h_plus[3], h_minus[3];
        for(int k=0; k<3; k++){
            h_plus[k]  = 0.5 * (UL[k] + phys_FL[k]/alpha);
            h_minus[k] = 0.5 * (UR[k] - phys_FR[k]/alpha);
        }

        pp_limiter(F_interface_plus, h_plus, F_interface_plus_re, eps1);
        pp_limiter(F_interface_minus, h_minus, F_interface_minus_re, eps2);

        for(int k=0; k<3; k++) {
            ws->Fhat[i*NVAR + k] = alpha * (F_interface_plus_re[k] - F_interface_minus_re[k]); 
        }
    }
}


// ================== 完备的精确黎曼求解器与有限体积初值积分 ==================

// 计算波函数 f_K 及其导数 df_K
static double eval_f_K(double p, double rho_K, double p_K, double c_K, double *df_K) {
    double f_val;
    if (p <= p_K) { // 稀疏波
        double pratio = p / p_K;
        double a = (GAMMA - 1.0) / (2.0 * GAMMA);
        f_val = (2.0 * c_K / (GAMMA - 1.0)) * (pow(pratio, a) - 1.0);
        *df_K = (1.0 / (rho_K * c_K)) * pow(pratio, -(GAMMA + 1.0) / (2.0 * GAMMA));
    } else {        // 激波
        double A = 2.0 / ((GAMMA + 1.0) * rho_K);
        double B = (GAMMA - 1.0) / (GAMMA + 1.0) * p_K;
        double sqrt_term = sqrt(A / (p + B));
        f_val = (p - p_K) * sqrt_term;
        *df_K = sqrt_term * (1.0 - (p - p_K) / (2.0 * (p + B)));
    }
    return f_val;
}

// 严格空间采样，返回精确物理量
static void exact_riemann_solution(double x, double t, double x_discontinuity, 
                                   double rho_L, double u_L, double p_L, 
                                   double rho_R, double u_R, double p_R,
                                   double *rho_out, double *u_out, double *p_out) {
    double c_L = sqrt(GAMMA * p_L / rho_L);
    double c_R = sqrt(GAMMA * p_R / rho_R);
    
    // TRRS (Two-Rarefaction Riemann Solver) 初值猜测
    double z = (GAMMA - 1.0) / (2.0 * GAMMA);
    double p_L_pow = pow(p_L, z);
    double p_R_pow = pow(p_R, z);
    double p_star = pow((c_L + c_R - 0.5 * (GAMMA - 1.0) * (u_R - u_L)) / (c_L / p_L_pow + c_R / p_R_pow), 1.0 / z);
    if (p_star < 1e-6 || isnan(p_star)) p_star = 0.5 * (p_L + p_R);

    // Newton-Raphson 代数方程求根
    double tol = 1e-7;
    for (int iter = 0; iter < 50; iter++) {
        double df_L, df_R;
        double f_val = eval_f_K(p_star, rho_L, p_L, c_L, &df_L) + eval_f_K(p_star, rho_R, p_R, c_R, &df_R) + (u_R - u_L);
        double df_val = df_L + df_R;
        
        double dp = -f_val / df_val;
        double p_new = p_star + dp;
        
        // 代数求根界限，非偏微分方程强制截断
        if (p_new < 1e-6) p_new = 1e-6; 
        
        double err = fabs(p_new - p_star) / (p_star + 1e-6);
        p_star = p_new;
        if (err < tol) break;
    }
    
    // 计算中间速度 u_star
    double df_dummy1, df_dummy2;
    double u_star = 0.5 * (u_L + u_R + eval_f_K(p_star, rho_R, p_R, c_R, &df_dummy1) - eval_f_K(p_star, rho_L, p_L, c_L, &df_dummy2));
    
    double S = (x - x_discontinuity) / t;
    
    // 全波系判定
    if (S <= u_star) { 
        if (p_star <= p_L) { // 左稀疏波
            double S_HL = u_L - c_L;
            if (S <= S_HL) { 
                *rho_out = rho_L; *u_out = u_L; *p_out = p_L; 
            } else {
                double c_star_L = c_L * pow(p_star / p_L, (GAMMA - 1.0) / (2.0 * GAMMA));
                double S_TL = u_star - c_star_L;
                if (S >= S_TL) {
                    *rho_out = rho_L * pow(p_star / p_L, 1.0 / GAMMA);
                    *u_out = u_star; *p_out = p_star;
                } else { // 稀疏波扇区内部
                    double u_fan = (2.0 / (GAMMA + 1.0)) * (c_L + 0.5 * (GAMMA - 1.0) * u_L + S);
                    double c_fan = (2.0 / (GAMMA + 1.0)) * (c_L + 0.5 * (GAMMA - 1.0) * (u_L - S));
                    *rho_out = rho_L * pow(c_fan / c_L, 2.0 / (GAMMA - 1.0));
                    *u_out = u_fan;
                    *p_out = p_L * pow(c_fan / c_L, 2.0 * GAMMA / (GAMMA - 1.0));
                }
            }
        } else { // 左激波
            double S_L = u_L - c_L * sqrt(((GAMMA + 1.0) / (2.0 * GAMMA)) * (p_star / p_L) + ((GAMMA - 1.0) / (2.0 * GAMMA)));
            if (S <= S_L) {
                *rho_out = rho_L; *u_out = u_L; *p_out = p_L;
            } else {
                double pratio = p_star / p_L;
                *rho_out = rho_L * (pratio + (GAMMA - 1.0) / (GAMMA + 1.0)) / (pratio * (GAMMA - 1.0) / (GAMMA + 1.0) + 1.0);
                *u_out = u_star; *p_out = p_star;
            }
        }
    } else { 
        if (p_star <= p_R) { // 右稀疏波
            double S_HR = u_R + c_R;
            if (S >= S_HR) {
                *rho_out = rho_R; *u_out = u_R; *p_out = p_R;
            } else {
                double c_star_R = c_R * pow(p_star / p_R, (GAMMA - 1.0) / (2.0 * GAMMA));
                double S_TR = u_star + c_star_R;
                if (S <= S_TR) {
                    *rho_out = rho_R * pow(p_star / p_R, 1.0 / GAMMA);
                    *u_out = u_star; *p_out = p_star;
                } else { // 稀疏波扇区内部
                    double u_fan = (2.0 / (GAMMA + 1.0)) * (-c_R + 0.5 * (GAMMA - 1.0) * u_R + S);
                    double c_fan = (2.0 / (GAMMA + 1.0)) * (c_R - 0.5 * (GAMMA - 1.0) * (u_R - S));
                    *rho_out = rho_R * pow(c_fan / c_R, 2.0 / (GAMMA - 1.0));
                    *u_out = u_fan;
                    *p_out = p_R * pow(c_fan / c_R, 2.0 * GAMMA / (GAMMA - 1.0));
                }
            }
        } else { // 右激波
            double S_R = u_R + c_R * sqrt(((GAMMA + 1.0) / (2.0 * GAMMA)) * (p_star / p_R) + ((GAMMA - 1.0) / (2.0 * GAMMA)));
            if (S >= S_R) {
                *rho_out = rho_R; *u_out = u_R; *p_out = p_R;
            } else {
                double pratio = p_star / p_R;
                *rho_out = rho_R * (pratio + (GAMMA - 1.0) / (GAMMA + 1.0)) / (pratio * (GAMMA - 1.0) / (GAMMA + 1.0) + 1.0);
                *u_out = u_star; *p_out = p_star;
            }
        }
    }
}

// 采用高斯积分在控制体积内求精确物理量的守恒均值
void init_condition_exact_riemann_fv(double *U, int Nx, double xL, double dx, double t_init) {
    const int N_GQ5 = 5;
    const double nodes_GQ5[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
    const double weights_GQ5[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

    double x_discontinuity = 0.3;
    double rho_L = 10000.0, u_L = 0.0, p_L = 10000.0;
    double rho_R = 1.0,     u_R = 0.0, p_R = 1.0;

    for (int i = 0; i < Nx; i++) {
        double xc = xL + (i + 0.5) * dx;
        
        double sum_rho = 0.0, sum_rhou = 0.0, sum_E = 0.0;
        
        // 在区间 [x_{i-1/2}, x_{i+1/2}] 进行积分
        for (int q = 0; q < N_GQ5; q++) {
            double xi = nodes_GQ5[q];
            double x_phys = xc + xi * (dx / 2.0);
            double w = weights_GQ5[q];
            
            double rho, u, p;
            exact_riemann_solution(x_phys, t_init, x_discontinuity, rho_L, u_L, p_L, rho_R, u_R, p_R, &rho, &u, &p);
            
            double rhou = rho * u;
            double E = p / (GAMMA - 1.0) + 0.5 * rho * u * u;
            
            sum_rho  += w * rho;
            sum_rhou += w * rhou;
            sum_E    += w * E;
        }
        
        // 将局部积分缩放到区间全长。
        // 原本需要乘以 (dx / 2.0) 再除以 dx，化简后常数系数为 0.5
        U[i * NVAR + 0] = 0.5 * sum_rho * lambdacoff;
        U[i * NVAR + 1] = 0.5 * sum_rhou * lambdacoff;
        U[i * NVAR + 2] = 0.5 * sum_E * lambdacoff;
    }
}

void residual_1D(double *U, double *R, SolverWorkspace *ws, double dx) {
    apply_BC_1D(U, ws);
    compute_Fhat_1D(ws);

    double inv_dx = 1.0 / dx;
    for (int i = 0; i < ws->Nx; i++) {
        for (int k = 0; k < NVAR; k++) {
            R[i * NVAR + k] = -(ws->Fhat[(i + 1) * NVAR + k] - ws->Fhat[i * NVAR + k]) * inv_dx;
        }
    }
}

void output_data(double *U, int Nx, double *xc) {
    FILE *fp = fopen("result_1D.txt", "w");
    if (!fp) return;
    fprintf(fp, "x rho u p\n");
    for (int i = 0; i < Nx; i++) {
        double rho = U[i*NVAR];
        double u = U[i*NVAR + 1] / rho;
        double p = pressure(&U[i*NVAR]);
    
        fprintf(fp, "%.12e %.12e %.12e %.12e\n", xc[i], rho, u, p);
    }
    fclose(fp);
    printf("Saved to result_1D.txt\n");
}

int main(void) {   
    _control87(0, _MCW_EM); 
    _control87(~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW), _MCW_EM);
    
    const double xL = 0, xR = 1;
    const double CFL = 0.1;
    const double t_final = 0.13;

    int Nx = 801;
    double dx = (xR - xL) / Nx;

    double *U  = (double*)malloc(Nx * NVAR * sizeof(double));
    double *U1 = (double*)malloc(Nx * NVAR * sizeof(double));
    double *U2 = (double*)malloc(Nx * NVAR * sizeof(double));
    double *L  = (double*)malloc(Nx * NVAR * sizeof(double));
    double *xc = (double*)malloc(Nx * sizeof(double));

    SolverWorkspace *ws = create_workspace(Nx);

    // Sod Shock Tube Initial Conditions
    for (int i = 0; i < Nx; ++i) {
        xc[i] = xL + (i + 0.5) * dx;
        double rho, u, E,p;
        if (fabs(xc[i] )< 0.3) { rho = 10000.0; u = 0; p =10000.0 ; } 
        else             { rho = 1.0; u = 0; p = 1.0; }
        E=p/ (GAMMA - 1.0) + 0.5 * (u * u) * rho;
        U[i*NVAR + 0] = rho*lambdacoff;
        U[i*NVAR + 1] = rho * u*lambdacoff;
        U[i*NVAR + 2] = E*lambdacoff;
    }

   // 【修改为】：提前一个微小时间步的有限体积积分初值
    double t = 0.01; 
    
    // 初始化坐标点数组供后续输出使用
    for (int i = 0; i < Nx; ++i) {
        xc[i] = xL + (i + 0.5) * dx;
    }

    // 调用精确解初值生成器填充控制体积均值
    init_condition_exact_riemann_fv(U, Nx, xL, dx, t);

    int it = 0;
    while (t < t_final) {
        double lambda_max = 0.0;

        for (int i = 0; i < Nx; i++) {
            double rho = U[i*NVAR], rhou = U[i*NVAR+1], E = U[i*NVAR+2];
            double u = rhou / rho;
            double p = pressure(&U[i*NVAR]);
            double c = sqrt(GAMMA * p / rho);
            double lam = fabs(u) + c;
            if (lam > lambda_max) lambda_max = lam;
        }

        double dt = CFL * dx / lambda_max;
        if (t + dt > t_final) dt = t_final - t;

        residual_1D(U, L, ws, dx);
        for (int i = 0; i < Nx*NVAR; ++i) U1[i] = U[i] + dt * L[i];

        residual_1D(U1, L, ws, dx);
        for (int i = 0; i < Nx*NVAR; ++i) 
            U2[i] = 0.75 * U[i] + 0.25 * (U1[i] + dt * L[i]);

        residual_1D(U2, L, ws, dx);
        for (int i = 0; i < Nx*NVAR; ++i) 
            U[i] = (1.0/3.0) * U[i] + (2.0/3.0) * (U2[i] + dt * L[i]);

        t += dt; it++;
      printf("it=%4d, t=%.5f, dt=%.3e\n", it, t, dt);
    }

    output_data(U, Nx, xc);
    for (int i = 0; i < Nx; ++i) {
        xc[i] = xL + (i + 0.5) * dx;
        double rho, u, E;
        if (fabs(xc[i]) < 1e-3) { printf("1!,%10f\n",xc[i]);rho = 1.0; u = 0; E = 3200000/dx; } 
        else             { rho = 1; u = 0; E = 1e-12; }
        
        U[i*NVAR + 0] = rho*lambdacoff;
        U[i*NVAR + 1] = rho * u*lambdacoff;
        U[i*NVAR + 2] = E*lambdacoff;
    }
    free_workspace(ws);
    free(U); free(U1); free(U2); free(L); free(xc);
    return 0;
}