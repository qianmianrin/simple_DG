#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define M_PI 3.14159265358979323846
/* ============================================================
 * 基本全局参数
 * ============================================================ */
#define Nx 400
#define Ny 400
#define xL 0.0
#define xR 1.0
#define yL 0.0
#define yR 1.0
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)
#define gamma 1.4
#define NUM_VARS 4 // 欧拉方程的变量组数：rho, rho u, rho v, E

/* 时间推进参数 */
#define CFL   0.1
#define T_END 0.3

/* ============================================================
 * 数据结构: 将标量扩展为长度为 4 的数组
 * ============================================================ */
typedef struct {
    double U[NUM_VARS][9]; 
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

double r_quad[9], s_quad[9], w_quad[9]; 

/* ============================================================
 * 预计算矩阵
 * ============================================================ */
#define N9 9
double M_inv[N9][N9];
double Sr[N9][N9];
double Ss[N9][N9];
double face_interp_T[3][N9];
double face_interp_B[3][N9];
double face_interp_L[3][N9];
double face_interp_R[3][N9];

/* ============================================================
 * 矩阵辅助函数
 * ============================================================ */
void mat_inv_9(double A[N9][N9], double Ainv[N9][N9]) {
    double temp[N9][2 * N9];
    for (int i = 0; i < N9; i++) {
        for (int j = 0; j < N9; j++) {
            temp[i][j] = A[i][j];
            temp[i][j + N9] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (int i = 0; i < N9; i++) {
        double p = temp[i][i];
        int row = i;
        for (int k = i + 1; k < N9; k++) 
            if (fabs(temp[k][i]) > fabs(p)) { p = temp[k][i]; row = k; }
        
        for (int j = 0; j < 2 * N9; j++) { double t = temp[i][j]; temp[i][j] = temp[row][j]; temp[row][j] = t; }
        
        double div = temp[i][i];
        for (int j = 0; j < 2 * N9; j++) temp[i][j] /= div;
        for (int k = 0; k < N9; k++) {
            if (k != i) {
                double factor = temp[k][i];
                for (int j = 0; j < 2 * N9; j++) temp[k][j] -= factor * temp[i][j];
            }
        }
    }
    for (int i = 0; i < N9; i++)
        for (int j = 0; j < N9; j++) Ainv[i][j] = temp[i][j + N9];
}

static inline double lagrange_1d(int i, double x) {
    double val = 1.0;
    for (int j = 0; j < 3; j++) {
        if (i == j) continue;
        val *= (x - nodes1D[j]) / (nodes1D[i] - nodes1D[j]);
    }
    return val;
}

static inline double d_lagrange_1d(int i, double x) {
    double deriv = 0.0;
    for (int j = 0; j < 3; j++) {
        if (i == j) continue;
        double term = 1.0 / (nodes1D[i] - nodes1D[j]);
        for (int k = 0; k < 3; k++) {
            if (k == i || k == j) continue;
            term *= (x - nodes1D[k]) / (nodes1D[i] - nodes1D[k]);
        }
        deriv += term;
    }
    return deriv;
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
    double M[N9][N9]  = {{0}};
    for (int q = 0; q < 9; q++) {
        double r = r_quad[q], s = s_quad[q], w = w_quad[q];
        for (int i = 0; i < N9; i++) {
            int ni = i % 3, mi = i / 3;
            double phi_i = lagrange_1d(ni, r) * lagrange_1d(mi, s);
            for (int j = 0; j < N9; j++) {
                int nj = j % 3, mj = j / 3;
                double phi_j      = lagrange_1d(nj, r)   * lagrange_1d(mj, s);
                double dphi_j_dr  = d_lagrange_1d(nj, r) * lagrange_1d(mj, s);
                double dphi_j_ds  = lagrange_1d(nj, r)   * d_lagrange_1d(mj, s);
                M[i][j]  += w * phi_i * phi_j;
                Sr[i][j] += w * phi_i * dphi_j_dr;
                Ss[i][j] += w * phi_i * dphi_j_ds;
            }
        }
    }
    for (int k = 0; k < 3; k++) {
        double nk = nodes1D[k];
        for (int j = 0; j < N9; j++) {
            int n = j % 3, m = j / 3;
            face_interp_T[k][j] = lagrange_1d(n, nk) * lagrange_1d(m,  1.0);
            face_interp_B[k][j] = lagrange_1d(n, nk) * lagrange_1d(m, -1.0);
            face_interp_L[k][j] = lagrange_1d(n, -1.0) * lagrange_1d(m, nk);
            face_interp_R[k][j] = lagrange_1d(n,  1.0) * lagrange_1d(m, nk);
        }
    }
    mat_inv_9(M, M_inv);
}

/* ============================================================
 * 欧拉方程物理通量与状态方程
 * ============================================================ */
static inline double calc_pressure(double rho, double rhou, double rhov, double E) {
    return (gamma - 1.0) * (E - 0.5 * (rhou * rhou + rhov * rhov) / rho);
}

static inline void euler_flux(double U[NUM_VARS], double F[NUM_VARS], double G[NUM_VARS]) {
    double rho = U[0];
    double rhou = U[1];
    double rhov = U[2];
    double E = U[3];
    
    double u = rhou / rho;
    double v = rhov / rho;
    double p = calc_pressure(rho, rhou, rhov, E);

    // x-direction flux F(U)
    F[0] = rhou;
    F[1] = rhou * u + p;
    F[2] = rhou * v;
    F[3] = u * (E + p);

    // y-direction flux G(U)
    G[0] = rhov;
    G[1] = rhou * v;
    G[2] = rhov * v + p;
    G[3] = v * (E + p);
}

static inline double max_speed() { 
    double max_v = 0.0;
    for (int i = 0; i < Nx * Ny; i++) {
        for (int q = 0; q < 9; q++) {
            double rho = Mesh[i].U[0][q];
            double rhou = Mesh[i].U[1][q];
            double rhov = Mesh[i].U[2][q];
            double E = Mesh[i].U[3][q];
            
            double u = rhou / rho;
            double v = rhov / rho;
            double p = calc_pressure(rho, rhou, rhov, E);
            double c = sqrt(gamma * p / rho);
            
            // 系统最大特征值发生在 |u|+c 或 |v|+c 处
            max_v = fmax(max_v, fabs(u) + c);
            max_v = fmax(max_v, fabs(v) + c);
        }
    }
    return (max_v < 1e-9) ? 1.0 : max_v; 
}

/* 向量化的 LF 数值通量 */
static inline void lf_flux_vector(double UL[NUM_VARS], double UR[NUM_VARS], double nx, double ny, double flux_res[NUM_VARS]) {
    double FL[NUM_VARS], GL[NUM_VARS];
    double FR[NUM_VARS], GR[NUM_VARS];
    
    euler_flux(UL, FL, GL);
    euler_flux(UR, FR, GR);

    // 计算左侧波速
    double pL = calc_pressure(UL[0], UL[1], UL[2], UL[3]);
    double cL = sqrt(gamma * pL / UL[0]);
    double vnL = (UL[1]/UL[0])*nx + (UL[2]/UL[0])*ny;
    double alphaL = fabs(vnL) + cL;

    // 计算右侧波速
    double pR = calc_pressure(UR[0], UR[1], UR[2], UR[3]);
    double cR = sqrt(gamma * pR / UR[0]);
    double vnR = (UR[1]/UR[0])*nx + (UR[2]/UR[0])*ny;
    double alphaR = fabs(vnR) + cR;

    double alpha = fmax(alphaL, alphaR);

    for(int v = 0; v < NUM_VARS; v++) {
        flux_res[v] = 0.5 * (FL[v]*nx + GL[v]*ny + FR[v]*nx + GR[v]*ny) - 0.5 * alpha * (UR[v] - UL[v]);
    }
}

/* ============================================================
 * 初始条件: 2D 黎曼问题 (四象限激波管)
 * ============================================================ */
void init_condition(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];
            double x_center = xL + (ii + 0.5) * dx;
            double y_center = yL + (jj + 0.5) * dy;

            for (int q = 0; q < 9; q++) {
                double x = x_center + 0.5 * dx * r_quad[q];
                double y = y_center + 0.5 * dy * s_quad[q];

                double rho, u, v, p;
                if (x >= 0.5 && y >= 0.5) { // Q1
                    rho = 1.5; u = 0.0; v = 0.0; p = 1.5;
                } else if (x < 0.5 && y >= 0.5) { // Q2(0.5323,1.206,0.0,0.3)
                    rho = 0.5323; u = 1.206; v = 0.0; p = 0.3;
                } else if (x < 0.5 && y < 0.5) { // Q3⎨

                    rho = 0.138; u = 1.206; v = 1.206; p =0.029;
                } else { // Q4(0.5323,0.0, 1.206,0.3)
                    rho = 0.5323; u = 0.0; v = 1.206; p = 0.3;
                }

                cell->U[0][q] = rho;
                cell->U[1][q] = rho * u;
                cell->U[2][q] = rho * v;
                cell->U[3][q] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v);
            }
        }
    }
}

/* ============================================================
 * 边界值插值（体节点 -> 面 Gauss 点）
 * ============================================================ */
void compute_face_values(Element *cell) {
    for (int var = 0; var < NUM_VARS; var++) {
        for (int k = 0; k < 3; k++) {
            double vt = 0.0, vb = 0.0, vl = 0.0, vr = 0.0;
            for (int j = 0; j < N9; j++) {
                double u = cell->U[var][j];
                vt += face_interp_T[k][j] * u;
                vb += face_interp_B[k][j] * u;
                vl += face_interp_L[k][j] * u;
                vr += face_interp_R[k][j] * u;
            }
            cell->face_top[var][k]    = vt;
            cell->face_bottom[var][k] = vb;
            cell->face_left[var][k]   = vl;
            cell->face_right[var][k]  = vr;
        }
    }
}

void boundary_value(void) {
    for (int idx = 0; idx < Nx * Ny; idx++)
        compute_face_values(&Mesh[idx]);
}

/* ============================================================
 * Ghost Cell 边界条件 (透射边界)
 * ============================================================ */
void apply_ghost_cells(void) {
    for (int i = 0; i < Nx; i++) {
        for (int var = 0; var < NUM_VARS; var++) {
            for (int k = 0; k < 3; k++) {
                Ghost_bottom[i].face_top[var][k] = Mesh[0 * Nx + i].face_bottom[var][k];
                Ghost_top[i].face_bottom[var][k] = Mesh[(Ny - 1) * Nx + i].face_top[var][k];
            }
        }
    }
    for (int j = 0; j < Ny; j++) {
        for (int var = 0; var < NUM_VARS; var++) {
            for (int k = 0; k < 3; k++) {
                Ghost_left[j].face_right[var][k] = Mesh[j * Nx + 0].face_left[var][k];
                Ghost_right[j].face_left[var][k] = Mesh[j * Nx + (Nx - 1)].face_right[var][k];
            }
        }
    }
}

static void get_neighbor_face(int ii, int jj, int face, double u_plus[NUM_VARS][3]) {
    for(int var = 0; var < NUM_VARS; var++) {
        if (face == 0) { 
            if (jj == 0) { for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_bottom[ii].face_top[var][k]; } 
            else { for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[(jj - 1) * Nx + ii].face_top[var][k]; }
        } else if (face == 1) { 
            if (jj == Ny - 1) { for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_top[ii].face_bottom[var][k]; } 
            else { for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[(jj + 1) * Nx + ii].face_bottom[var][k]; }
        } else if (face == 2) { 
            if (ii == 0) { for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_left[jj].face_right[var][k]; } 
            else { for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[jj * Nx + (ii - 1)].face_right[var][k]; }
        } else { 
            if (ii == Nx - 1) { for (int k = 0; k < 3; k++) u_plus[var][k] = Ghost_right[jj].face_left[var][k]; } 
            else { for (int k = 0; k < 3; k++) u_plus[var][k] = Mesh[jj * Nx + (ii + 1)].face_left[var][k]; }
        }
    }
}

/* ============================================================
 * 限制器 (分量限制法)
 * ============================================================ */
static void cell_avg(int idx, double Ubar[NUM_VARS])
{
    double sw = 0.0;
    int var;
    for (var = 0; var < NUM_VARS; var++) Ubar[var] = 0.0;
    for (int q = 0; q < 9; q++) {
        double w = w_quad[q];
        sw += w;
        for (var = 0; var < NUM_VARS; var++)
            Ubar[var] += w * Mesh[idx].U[var][q];
    }
    for (var = 0; var < NUM_VARS; var++) Ubar[var] /= sw;
}

static double minmod_tvb(double a, double b, double c, double h)
{
    static const double M_TVB = 50.0;
    if (fabs(a) <= M_TVB * h * h) return a;
    if (a > 0.0 && b > 0.0 && c > 0.0) return fmin(a, fmin(b, c));
    if (a < 0.0 && b < 0.0 && c < 0.0) return fmax(a, fmax(b, c));
    return 0.0;
}


static void build_eigen_x(double rho, double u, double v, double p,
                           double R[NUM_VARS][NUM_VARS],
                           double L[NUM_VARS][NUM_VARS])
{
    double c   = sqrt(gamma * p / rho);
    double H   = (p/(gamma-1.0) + 0.5*rho*(u*u+v*v) + p) / rho;
    double q2  = 0.5*(u*u + v*v);
    double c2  = c*c;
    double gm1 = gamma - 1.0;
    double a1  = gm1 / (2.0*c2);
    double b   = gm1 / c2;

    /* R_x: 列为右特征向量 */
    R[0][0]=1.0;  R[1][0]=u-c;  R[2][0]=v;   R[3][0]=H-u*c;
    R[0][1]=1.0;  R[1][1]=u;    R[2][1]=v;   R[3][1]=q2;
    R[0][2]=0.0;  R[1][2]=0.0;  R[2][2]=1.0; R[3][2]=v;
    R[0][3]=1.0;  R[1][3]=u+c;  R[2][3]=v;   R[3][3]=H+u*c;

    /* L_x = R_x^{-1}: 行为左特征向量 */
    L[0][0]= a1*q2 + u/(2.0*c); L[0][1]=-a1*u - 1.0/(2.0*c); L[0][2]=-a1*v; L[0][3]=a1;
    L[1][0]= 1.0 - b*q2;        L[1][1]= b*u;                 L[1][2]= b*v;  L[1][3]=-b;
    L[2][0]=-v;                  L[2][1]= 0.0;                 L[2][2]= 1.0;  L[2][3]=0.0;
    L[3][0]= a1*q2 - u/(2.0*c); L[3][1]=-a1*u + 1.0/(2.0*c); L[3][2]=-a1*v; L[3][3]=a1;
}

static void build_eigen_y(double rho, double u, double v, double p,
                           double R[NUM_VARS][NUM_VARS],
                           double L[NUM_VARS][NUM_VARS])
{
    double c   = sqrt(gamma * p / rho);
    double H   = (p/(gamma-1.0) + 0.5*rho*(u*u+v*v) + p) / rho;
    double q2  = 0.5*(u*u + v*v);
    double c2  = c*c;
    double gm1 = gamma - 1.0;
    double a1  = gm1 / (2.0*c2);
    double b   = gm1 / c2;

    /* R_y: 列为右特征向量 (u <-> v 对称) */
    R[0][0]=1.0;  R[1][0]=u;   R[2][0]=v-c;  R[3][0]=H-v*c;
    R[0][1]=1.0;  R[1][1]=u;   R[2][1]=v;    R[3][1]=q2;
    R[0][2]=0.0;  R[1][2]=1.0; R[2][2]=0.0;  R[3][2]=u;
    R[0][3]=1.0;  R[1][3]=u;   R[2][3]=v+c;  R[3][3]=H+v*c;

    /* L_y = R_y^{-1} */
    L[0][0]= a1*q2 + v/(2.0*c); L[0][1]=-a1*u; L[0][2]=-a1*v - 1.0/(2.0*c); L[0][3]=a1;
    L[1][0]= 1.0 - b*q2;        L[1][1]= b*u;  L[1][2]= b*v;                 L[1][3]=-b;
    L[2][0]=-u;                  L[2][1]= 1.0;  L[2][2]= 0.0;                 L[2][3]=0.0;
    L[3][0]= a1*q2 - v/(2.0*c); L[3][1]=-a1*u; L[3][2]=-a1*v + 1.0/(2.0*c); L[3][3]=a1;
}

void apply_limiter(void)
{
    int ii, jj, var, m, n, q;

    for (jj = 0; jj < Ny; jj++) {
        for (ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;

            double Ubar[NUM_VARS];
            cell_avg(idx, Ubar);

            double UbarL[NUM_VARS], UbarR[NUM_VARS];
            double UbarB[NUM_VARS], UbarT[NUM_VARS];

            //不太对。为了更统一的处理边界，改为Gosh Cell
            if (ii > 0)    cell_avg(jj*Nx + ii-1,   UbarL);
            else           for(var=0;var<NUM_VARS;var++) UbarL[var]=Ubar[var];
            if (ii < Nx-1) cell_avg(jj*Nx + ii+1,   UbarR);
            else           for(var=0;var<NUM_VARS;var++) UbarR[var]=Ubar[var];
            if (jj > 0)    cell_avg((jj-1)*Nx + ii, UbarB);
            else           for(var=0;var<NUM_VARS;var++) UbarB[var]=Ubar[var];
            if (jj < Ny-1) cell_avg((jj+1)*Nx + ii, UbarT);
            else           for(var=0;var<NUM_VARS;var++) UbarT[var]=Ubar[var];

            double rho = Ubar[0];
            double u   = Ubar[1] / rho;
            double v   = Ubar[2] / rho;
            double p   = (gamma-1.0)*(Ubar[3]
                         - 0.5*(Ubar[1]*Ubar[1]+Ubar[2]*Ubar[2])/rho);
            if (rho <= 0.0 || p <= 0.0) continue;

            /* 两套特征矩阵 */
            double Rx[NUM_VARS][NUM_VARS], Lx[NUM_VARS][NUM_VARS];
            double Ry[NUM_VARS][NUM_VARS], Ly[NUM_VARS][NUM_VARS];
            build_eigen_x(rho, u, v, p, Rx, Lx);
            build_eigen_y(rho, u, v, p, Ry, Ly);

            /* 面中值 - 均值 (守恒空间斜率) */
            double dUx[NUM_VARS], dUy[NUM_VARS];
            for (var = 0; var < NUM_VARS; var++) {
                double uR_face = 0.0, uT_face = 0.0;
                for (int k = 0; k < N9; k++) {
                    uR_face += face_interp_R[1][k] * Mesh[idx].U[var][k];
                    uT_face += face_interp_T[1][k] * Mesh[idx].U[var][k];
                }
                dUx[var] = uR_face - Ubar[var];
                dUy[var] = uT_face - Ubar[var];
            }

            double dUx_R[NUM_VARS], dUx_L[NUM_VARS];
            double dUy_T[NUM_VARS], dUy_B[NUM_VARS];
            for (var = 0; var < NUM_VARS; var++) {
                dUx_R[var] = UbarR[var] - Ubar[var];
                dUx_L[var] = Ubar[var]  - UbarL[var];
                dUy_T[var] = UbarT[var] - Ubar[var];
                dUy_B[var] = Ubar[var]  - UbarB[var];
            }

            /* x 方向: 用 Lx 投影 */
            double dWx[NUM_VARS], dWx_R[NUM_VARS], dWx_L[NUM_VARS];
            for (m = 0; m < NUM_VARS; m++) {
                dWx[m]=0.0; dWx_R[m]=0.0; dWx_L[m]=0.0;
                for (n = 0; n < NUM_VARS; n++) {
                    dWx[m]   += Lx[m][n] * dUx[n];
                    dWx_R[m] += Lx[m][n] * dUx_R[n];
                    dWx_L[m] += Lx[m][n] * dUx_L[n];
                }
            }

            /* y 方向: 用 Ly 投影 */
            double dWy[NUM_VARS], dWy_T[NUM_VARS], dWy_B[NUM_VARS];
            for (m = 0; m < NUM_VARS; m++) {
                dWy[m]=0.0; dWy_T[m]=0.0; dWy_B[m]=0.0;
                for (n = 0; n < NUM_VARS; n++) {
                    dWy[m]   += Ly[m][n] * dUy[n];
                    dWy_T[m] += Ly[m][n] * dUy_T[n];
                    dWy_B[m] += Ly[m][n] * dUy_B[n];
                }
            }

            /* TVB minmod (各自特征空间) */
            double dWx_lim[NUM_VARS], dWy_lim[NUM_VARS];
            int need_limit = 0;
            for (m = 0; m < NUM_VARS; m++) {
                dWx_lim[m] = minmod_tvb(dWx[m], dWx_R[m], dWx_L[m], dx);
                dWy_lim[m] = minmod_tvb(dWy[m], dWy_T[m], dWy_B[m], dy);
                if (fabs(dWx[m]-dWx_lim[m]) > 1e-14 ||
                    fabs(dWy[m]-dWy_lim[m]) > 1e-14)
                    need_limit = 1;
            }
            if (!need_limit) continue;

            /* 变换回守恒空间: x用Rx, y用Ry */
            double dUx_lim[NUM_VARS], dUy_lim[NUM_VARS];
            for (m = 0; m < NUM_VARS; m++) {
                dUx_lim[m]=0.0; dUy_lim[m]=0.0;
                for (n = 0; n < NUM_VARS; n++) {
                    dUx_lim[m] += Rx[m][n] * dWx_lim[n];
                    dUy_lim[m] += Ry[m][n] * dWy_lim[n];
                }
            }

            /* 线性重构 */
            for (q = 0; q < 9; q++) {
                double r = r_quad[q], s = s_quad[q];
                for (var = 0; var < NUM_VARS; var++) {
                    Mesh[idx].U[var][q] = Ubar[var]
                                        + dUx_lim[var] * r
                                        + dUy_lim[var] * s;
                }
            }
        }
    }
}

/* ============================================================
 * 强形式 DG 右端项
 * ============================================================ */
void compute_rhs(double RHS[Nx * Ny][NUM_VARS][N9]) {
    boundary_value();
    apply_ghost_cells();

    double J = (dx * dy) / 4.0; 
    double Jx = dy / 2.0; 
    double Jy = dx / 2.0; 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double surf_res[NUM_VARS][N9] = {{0}}; 
            double u_plus[NUM_VARS][3];

            /* --- 步骤 A: 体积分贡献 --- */
            for (int i = 0; i < N9; i++) {
                for(int var = 0; var < NUM_VARS; var++) RHS[idx][var][i] = 0.0;
                
                double dfdr[NUM_VARS] = {0}, dgds[NUM_VARS] = {0};
                for (int j = 0; j < N9; j++) {
                    double F_val[NUM_VARS], G_val[NUM_VARS], U_j[NUM_VARS];
                    for(int v=0; v<NUM_VARS; v++) U_j[v] = cell->U[v][j];
                    euler_flux(U_j, F_val, G_val);
                    
                    for(int v=0; v<NUM_VARS; v++) {
                        dfdr[v] += Sr[i][j] * F_val[v];
                        dgds[v] += Ss[i][j] * G_val[v];
                    }
                }
                for(int v=0; v<NUM_VARS; v++) {
                    RHS[idx][v][i] = -(Jx * dfdr[v] + Jy * dgds[v]);
                }
            }

            /* --- 步骤 B: 边界积分贡献 --- */
            double U_minus[NUM_VARS], U_p[NUM_VARS], num_f[NUM_VARS], F_int[NUM_VARS], G_int[NUM_VARS];

            // 1. Bottom (s = -1, n = [0, -1])
            get_neighbor_face(ii, jj, 0, u_plus);
            for (int k = 0; k < 3; k++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_bottom[v][k]; U_p[v] = u_plus[v][k]; }
                euler_flux(U_minus, F_int, G_int);
                lf_flux_vector(U_minus, U_p, 0, -1, num_f);
                
                for(int v=0; v<NUM_VARS; v++) {
                    double diff = -G_int[v] - num_f[v]; 
                    for (int i = 0; i < N9; i++) 
                        surf_res[v][i] += (dx / 2.0) * weights1D[k] * diff * face_interp_B[k][i];
                }
            }

            // 2. Top (s = 1, n = [0, 1])
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int k = 0; k < 3; k++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_top[v][k]; U_p[v] = u_plus[v][k]; }
                euler_flux(U_minus, F_int, G_int);
                lf_flux_vector(U_minus, U_p, 0, 1, num_f);
                
                for(int v=0; v<NUM_VARS; v++) {
                    double diff = G_int[v] - num_f[v];
                    for (int i = 0; i < N9; i++) 
                        surf_res[v][i] += (dx / 2.0) * weights1D[k] * diff * face_interp_T[k][i];
                }
            }

            // 3. Left (r = -1, n = [-1, 0])
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int k = 0; k < 3; k++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_left[v][k]; U_p[v] = u_plus[v][k]; }
                euler_flux(U_minus, F_int, G_int);
                lf_flux_vector(U_minus, U_p, -1, 0, num_f);
                
                for(int v=0; v<NUM_VARS; v++) {
                    double diff = -F_int[v] - num_f[v];
                    for (int i = 0; i < N9; i++) 
                        surf_res[v][i] += (dy / 2.0) * weights1D[k] * diff * face_interp_L[k][i];
                }
            }

            // 4. Right (r = 1, n = [1, 0])
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int k = 0; k < 3; k++) {
                for(int v=0; v<NUM_VARS; v++) { U_minus[v] = cell->face_right[v][k]; U_p[v] = u_plus[v][k]; }
                euler_flux(U_minus, F_int, G_int);
                lf_flux_vector(U_minus, U_p, 1, 0, num_f);
                
                for(int v=0; v<NUM_VARS; v++) {
                    double diff = F_int[v] - num_f[v];
                    for (int i = 0; i < N9; i++) 
                        surf_res[v][i] += (dy / 2.0) * weights1D[k] * diff * face_interp_R[k][i];
                }
            }

            /* --- 步骤 C: 乘以质量矩阵之逆 --- */
            for(int v=0; v<NUM_VARS; v++) {
                for (int i = 0; i < N9; i++) RHS[idx][v][i] += surf_res[v][i];

                double RHS_temp[N9] = {0};
                for (int i = 0; i < N9; i++) {
                    for(int col = 0; col < N9; col++) {
                        RHS_temp[i] += (M_inv[i][col] * RHS[idx][v][col]);
                    }
                }
                for (int i = 0; i < N9; i++) RHS[idx][v][i] = RHS_temp[i] / J;
            }
        }
    }
}

/* ============================================================
 * 时间步长
 * ============================================================ */
static double compute_dt(void) {
    double speed = max_speed();
    double h     = (dx < dy) ? dx : dy;
    return CFL * h / speed;
}

/* ============================================================
 * SSP-RK3
 * ============================================================ */
static double U0[Nx * Ny][NUM_VARS][N9]; 
static double U1[Nx * Ny][NUM_VARS][N9]; 
static double U2[Nx * Ny][NUM_VARS][N9]; 
static double RHS_buf[Nx * Ny][NUM_VARS][N9];

void rk3_step(double dt, int nit) {
    int total = Nx * Ny;
    for (int i = 0; i < total; i++) memcpy(U0[i], Mesh[i].U, NUM_VARS * N9 * sizeof(double));

    /* Stage 1 */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int v = 0; v < NUM_VARS; v++)
            for (int q = 0; q < N9; q++)
                U1[i][v][q] = U0[i][v][q] + dt * RHS_buf[i][v][q];
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U1[i], NUM_VARS * N9 * sizeof(double));
    apply_limiter();

    /* Stage 2 */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int v = 0; v < NUM_VARS; v++)
            for (int q = 0; q < N9; q++)
                U2[i][v][q] = 0.75 * U0[i][v][q] + 0.25 * (U1[i][v][q] + dt * RHS_buf[i][v][q]);
    for (int i = 0; i < total; i++) memcpy(Mesh[i].U, U2[i], NUM_VARS * N9 * sizeof(double));
    apply_limiter();

    /* Stage 3 */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int v = 0; v < NUM_VARS; v++)
            for (int q = 0; q < N9; q++)
                Mesh[i].U[v][q] = (1.0 / 3.0) * U0[i][v][q] + (2.0 / 3.0) * (U2[i][v][q] + dt * RHS_buf[i][v][q]);
    apply_limiter();
}

void output_results(double t) {
    char filename[50];
    sprintf(filename, "result.dat");
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) return;
    
    // 导出原始变量: 密度, x速度, y速度, 压力
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Rho\", \"U\", \"V\", \"P\"\n");

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double xc = (ii + 0.5) * dx;
            double yc = (jj + 0.5) * dy;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    int n_idx = i * 3 + j;
                    double x_phys = xc + (dx / 2.0) * nodes1D[j];
                    double y_phys = yc + (dy / 2.0) * nodes1D[i];
                    
                    double rho = cell->U[0][n_idx];
                    double u   = cell->U[1][n_idx] / rho;
                    double v   = cell->U[2][n_idx] / rho;
                    double E   = cell->U[3][n_idx];
                    double p   = calc_pressure(rho, cell->U[1][n_idx], cell->U[2][n_idx], E);

                    fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", x_phys, y_phys, rho, u, v, p);
                }
            }
        }
    }
    fclose(fp);
    printf("Saved results to %s\n", filename);
}

int main(void) {
    init_quadrature();
    precompute_matrices();
    init_condition();

    double t  = 0.0;
    int   nit = 0;

    printf("%-10s  %-14s\n", "Step", "Time");
    printf("--------------------------\n");

    while (t < T_END) {
        double dt = compute_dt();
        if (t + dt > T_END) dt = T_END - t;

        rk3_step(dt, nit);
        t  += dt;
        nit++;

            printf("%-10d  %-14.6e\n", nit, t);
      
    }
    output_results(T_END);
    return 0;
}