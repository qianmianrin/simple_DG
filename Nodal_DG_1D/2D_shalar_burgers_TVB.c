#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define M_PI 3.14159265358979323846
/* ============================================================
 * 基本全局参数
 * ============================================================ */
#define Nx 30
#define Ny 30
#define xL 0.0
#define xR 1.0
#define yL 0.0
#define yR 1.0
#define dx ((xR - xL) / Nx)
#define dy ((yR - yL) / Ny)

/* 时间推进参数 */
#define CFL   0.1
#define T_END 0.1

/* ============================================================
 * 数据结构
 * ============================================================
*/
typedef struct {
    double U[9];
    double face_top[3];
    double face_bottom[3];
    double face_left[3];
    double face_right[3];
} Element;

Element Mesh[Nx * Ny]; /* Mesh[j*Nx + i] 对应第 (i,j) 单元 */

/* Ghost Cell：四条外边界各 Nx 或 Ny 个虚单元 */
Element Ghost_bottom[Nx]; /* j=-1  行 */
Element Ghost_top[Nx];    /* j=Ny  行 */
Element Ghost_left[Ny];   /* i=-1  列 */
Element Ghost_right[Ny];  /* i=Nx  列 */

/* ============================================================
 * 3 点 Gauss-Legendre 节点与权重（参考区间 [-1,1]）
 * ============================================================ */
static const double nodes1D[3]   = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = { 0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

double r_quad[9], s_quad[9], w_quad[9]; /* 9 点 2D 张量积 */

/* ============================================================
 * 预计算矩阵
 * ============================================================
*/
double M_inv[9][9];
double Sr[9][9];
double Ss[9][9];
double face_interp_T[3][9];
double face_interp_B[3][9];
double face_interp_L[3][9];
double face_interp_R[3][9];


/* ============================================================
 * 物理通量函数
 * ============================================================ */

/* x 方向通量 f(u) */
static inline double flux_f(double u) {
    return 0.5*u* u;
}

/* y 方向通量 g(u) */
static inline double flux_g(double u) {
    return 0.5*u * u;
}

/* ============================================================
 * 最大波速（用于 LF 通量稳定性参数 alpha）
 * ============================================================ */
static inline double max_speed() { double max_v = 0.0;
    for (int i = 0; i < Nx * Ny; i++) {
        for (int q = 0; q < 9; q++) {
            max_v = fmax(max_v, fabs(Mesh[i].U[q]));
        }
    }
    return (max_v < 1e-9) ? 1.0 : max_v; // 防止初值为0导致 dt 为无穷大
    }


/* ============================================================
 * LF（Lax-Friedrichs）数值通量
 * ============================================================ */
static inline double lf_flux(double u_minus, double u_plus,double nx,double ny) {
 double alpha = fmax(fabs(u_minus), fabs(u_plus));
    return 0.5 * (flux_f(u_minus) *nx+ flux_f(u_plus)*nx+flux_g(u_minus)*ny+flux_g(u_plus)*ny)
           - 0.5 * alpha * (u_plus - u_minus);
}



/* ============================================================
矩阵辅助函数
 * ============================================================ */

#define N9 9
/* 矩阵乘法：C = A * B */
void mat_mul_9(double A[N9][N9], double B[N9][N9], double C[N9][N9]) {
    for (int i = 0; i < N9; i++) {
        for (int j = 0; j < N9; j++) {
            C[i][j] = 0;
            for (int k = 0; k < N9; k++) C[i][j] += A[i][k] * B[k][j];
        }
    }
}

/* 矩阵求逆：高斯-约当消元法 */
void mat_inv_9(double A[N9][N9], double Ainv[N9][N9]) {
    double temp[N9][2 * N9];
    // 构造增广矩阵 [A | I]
    for (int i = 0; i < N9; i++) {
        for (int j = 0; j < N9; j++) {
            temp[i][j] = A[i][j];
            temp[i][j + N9] = (i == j) ? 1.0 : 0.0;
        }
    }
    // 消元
    for (int i = 0; i < N9; i++) {
        double p = temp[i][i];
        int row = i;
        for (int k = i + 1; k < N9; k++) // 选主元保证数值稳定性
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
    // 提取右侧逆矩阵
    for (int i = 0; i < N9; i++)
        for (int j = 0; j < N9; j++) Ainv[i][j] = temp[i][j + N9];
}

/* ============================================================
 * 1D Lagrange 基函数及其导数
 * ============================================================ */
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

/* ============================================================
 * 初始化 2D 张量积 Gauss 积分点
 * ============================================================ */
void init_quadrature(void) {
    for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) {
            int idx     = j * 3 + i;
            r_quad[idx] = nodes1D[i];
            s_quad[idx] = nodes1D[j];
            w_quad[idx] = weights1D[i] * weights1D[j];
        }
}

/* ============================================================
 * 预计算所有算子矩阵
 * ============================================================ */
void precompute_matrices(void) {
    double M[N9][N9]  = {{0}};

    double tmp[N9][N9];

    /* 体积分 */
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

    /* 面插值矩阵 */
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

    /* M^{-1} */
    mat_inv_9(M, M_inv);

  
}

/* ============================================================
 * 初始条件  u_0(x,y) = sin(2*pi*x) * sin(2*pi*y)
 * ============================================================ */
void init_condition(void) {
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];
            double x_center = xL + (ii + 0.5) * dx;
            double y_center = yL + (jj + 0.5) * dy;

            for (int q = 0; q < 9; q++) {
                // 获取每个 Gauss 积分点的物理坐标
                double x = x_center + 0.5 * dx * r_quad[q];
                double y = y_center + 0.5 * dy * s_quad[q];

                if (x > 0.5 && y > 0.5)      cell->U[q] = 0.5;  // Q1
                else if (x <= 0.5 && y > 0.5) cell->U[q] = -0.2; // Q2
                else if (x <= 0.5 && y <= 0.5) cell->U[q] = -1.0; // Q3
                else                          cell->U[q] = 0.8;  // Q4
            }
        }
    }
}

/* ============================================================
 * 边界值插值（体节点 -> 面 Gauss 点）
 * ============================================================ */
void compute_face_values(Element *cell) {
    for (int k = 0; k < 3; k++) {
        double vt = 0.0, vb = 0.0, vl = 0.0, vr = 0.0;
        for (int j = 0; j < N9; j++) {
            double u = cell->U[j];
            vt += face_interp_T[k][j] * u;
            vb += face_interp_B[k][j] * u;
            vl += face_interp_L[k][j] * u;
            vr += face_interp_R[k][j] * u;
        }
        cell->face_top[k]    = vt;
        cell->face_bottom[k] = vb;
        cell->face_left[k]   = vl;
        cell->face_right[k]  = vr;
    }
}

void boundary_value(void) {
    for (int idx = 0; idx < Nx * Ny; idx++)
        compute_face_values(&Mesh[idx]);
}

/* ============================================================
 * Ghost Cell 边界条件
 * ============================================================ */

void apply_ghost_cells(void) {
    /* --- 1. 底部边界 (j=0) 与 顶部边界 (j=Ny-1) --- */
    for (int i = 0; i < Nx; i++) {
        for (int k = 0; k < 3; k++) {
            // 底部 Ghost (j=-1): 
            // 它的“上表面”应该接到 Mesh[0][i] 的“下表面”
            // 流出条件下，外侧值直接等于内侧值
            Ghost_bottom[i].face_top[k] = Mesh[0 * Nx + i].face_bottom[k];

            // 顶部 Ghost (j=Ny): 
            // 它的“下表面”接到 Mesh[Ny-1][i] 的“上表面”
            Ghost_top[i].face_bottom[k] = Mesh[(Ny - 1) * Nx + i].face_top[k];
        }
    }

    /* --- 2. 左侧边界 (i=0) 与 右侧边界 (i=Nx-1) --- */
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < 3; k++) {
            // 左侧 Ghost (i=-1): 
            // 它的“右表面”接到 Mesh[j][0] 的“左表面”
            Ghost_left[j].face_right[k] = Mesh[j * Nx + 0].face_left[k];

            // 右侧 Ghost (i=Nx): 
            // 它的“左表面”接到 Mesh[j][Nx-1] 的“右表面”
            Ghost_right[j].face_left[k] = Mesh[j * Nx + (Nx - 1)].face_right[k];
        }
    }
}
/* ============================================================
 * 获取相邻单元在共享面上的值（u^+，即外侧值）
 * ============================================================ */
static void get_neighbor_face(int ii, int jj, int face,
                               double u_plus[3]) {
    int k;
    if (face == 0) { /* bottom: 邻居在 j-1 */
        if (jj == 0) {
            for (k = 0; k < 3; k++) u_plus[k] = Ghost_bottom[ii].face_top[k];
        } else {
            for (k = 0; k < 3; k++) u_plus[k] = Mesh[(jj - 1) * Nx + ii].face_top[k];
        }
    } else if (face == 1) { /* top: 邻居在 j+1 */
        if (jj == Ny - 1) {
            for (k = 0; k < 3; k++) u_plus[k] = Ghost_top[ii].face_bottom[k];
        } else {
            for (k = 0; k < 3; k++) u_plus[k] = Mesh[(jj + 1) * Nx + ii].face_bottom[k];
        }
    } else if (face == 2) { /* left: 邻居在 i-1 */
        if (ii == 0) {
            for (k = 0; k < 3; k++) u_plus[k] = Ghost_left[jj].face_right[k];
        } else {
            for (k = 0; k < 3; k++) u_plus[k] = Mesh[jj * Nx + (ii - 1)].face_right[k];
        }
    } else { /* right: 邻居在 i+1 */
        if (ii == Nx - 1) {
            for (k = 0; k < 3; k++) u_plus[k] = Ghost_right[jj].face_left[k];
        } else {
            for (k = 0; k < 3; k++) u_plus[k] = Mesh[jj * Nx + (ii + 1)].face_left[k];
        }
    }
}

/* ============================================================
 * 限制器
 * ============================================================ */
double get_cell_avg(int idx) {
    double sum = 0;
    for (int q = 0; q < 9; q++) sum += w_quad[q] * Mesh[idx].U[q];
    return sum / 4.0; // 参考区间面积为 4
}
double minmod(double a, double b, double c,double delta) {
    double M = 20.0; // TVB 参数，根据激波强度调整
    if (fabs(a) <= M * delta * delta) return a;
    if (a > 0 && b > 0 && c > 0) return fmin(a, fmin(b, c));
    if (a < 0 && b < 0 && c < 0) return fmax(a, fmax(b, c));
    return 0.0;
}

void apply_limiter() {
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            int idx = j * Nx + i;
            double u_bar = get_cell_avg(idx);
            
            // 获取相邻单元平均值 (处理边界逻辑同 get_neighbor_face)
            double u_bar_L = (i > 0) ? get_cell_avg(j * Nx + i - 1) : u_bar; 
            double u_bar_R = (i < Nx - 1) ? get_cell_avg(j * Nx + i + 1) : u_bar;
            double u_bar_B = (j > 0) ? get_cell_avg((j - 1) * Nx + i) : u_bar;
            double u_bar_T = (j < Ny - 1) ? get_cell_avg((j + 1) * Nx + i) : u_bar;

            // 1. X 方向限制
            // 计算右边界中点值与均值的差
            double u_R = 0; 
            for(int k=0; k<9; k++) u_R += face_interp_R[1][k] * Mesh[idx].U[k];
            double del_x = u_R - u_bar;
            double del_x_mod = minmod(del_x, u_bar_R - u_bar, u_bar - u_bar_L,dx);

            // 2. Y 方向限制
            double u_T = 0;
            for(int k=0; k<9; k++) u_T += face_interp_T[1][k] * Mesh[idx].U[k];
            double del_y = u_T - u_bar;
            double del_y_mod = minmod(del_y, u_bar_T - u_bar, u_bar - u_bar_B,dy);

            // 3. 如果发生了限制，则将多项式重构为线性：
            // u(r,s) = u_bar + del_x_mod * r + del_y_mod * s
            if (fabs(del_x - del_x_mod) > 1e-10 || fabs(del_y - del_y_mod) > 1e-10) {
                for (int q = 0; q < 9; q++) {
                    Mesh[idx].U[q] = u_bar + del_x_mod * r_quad[q] + del_y_mod * s_quad[q];
                }
            }
        }
    }
}

/* ============================================================
 * 强形式 DG 右端项
 * ============================================================ */

void compute_rhs(double RHS[Nx * Ny][N9]) {
    /* 1. 预处理：计算面值并填充 Ghost Cell */
    boundary_value();
    apply_ghost_cells();

    // 几何因子预计算
    double J = (dx * dy) / 4.0; 
    // 体积分所需的链式法则系数缩放，此处直接给出化简结果
    // - \int (\partial f / \partial x) \phi_i dxdy = - (\Delta y / 2) * \int (\partial f / \partial r) \phi_i drds
    double Jx = dy / 2.0; 
    double Jy = dx / 2.0; 

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];
            double surf_res[N9] = {0}; 
            double u_plus[3];

            /* --- 步骤 A: 体积分贡献 (Volume Term) --- */
            for (int i = 0; i < N9; i++) {
                double dfdr = 0.0, dgds = 0.0;
                for (int j = 0; j < N9; j++) {
                    dfdr += Sr[i][j] * flux_f(cell->U[j]);
                    dgds += Ss[i][j] * flux_g(cell->U[j]);
                }
                // 体积分投影
                RHS[idx][i] = -(Jx * dfdr + Jy * dgds);
            }

            /* --- 步骤 B: 边界积分贡献 (Lifting Term) --- */
            // 1. Bottom (s = -1, n = [0, -1])
            get_neighbor_face(ii, jj, 0, u_plus);
            for (int k = 0; k < 3; k++) {
                double uM = cell->face_bottom[k];
                double uP = u_plus[k];
                double fn_int = -flux_g(uM); 
                // 严格的 LF 通量公式：0.5*(F_L*n + F_R*n) - 0.5*alpha*(u_R - u_L)
                double fn_num = lf_flux(uM,uP,0,-1);
                double diff = fn_int - fn_num;
                for (int i = 0; i < N9; i++) 
                    surf_res[i] += (dx / 2.0) * weights1D[k] * diff * face_interp_B[k][i];
            }

            // 2. Top (s = 1, n = [0, 1])
            get_neighbor_face(ii, jj, 1, u_plus);
            for (int k = 0; k < 3; k++) {
                double uM = cell->face_top[k];
                double uP = u_plus[k];
                double fn_int = flux_g(uM);
                double fn_num = lf_flux(uM,uP,0,1);
                double diff = fn_int - fn_num;
                for (int i = 0; i < N9; i++) 
                    surf_res[i] += (dx / 2.0) * weights1D[k] * diff * face_interp_T[k][i];
            }

            // 3. Left (r = -1, n = [-1, 0])
            get_neighbor_face(ii, jj, 2, u_plus);
            for (int k = 0; k < 3; k++) {
                double uM = cell->face_left[k];
                double uP = u_plus[k];
                double fn_int = -flux_f(uM);
                double fn_num = lf_flux(uM,uP,-1,0);
                double diff = fn_int - fn_num;
                for (int i = 0; i < N9; i++) 
                    surf_res[i] += (dy / 2.0) * weights1D[k] * diff * face_interp_L[k][i];
            }

            // 4. Right (r = 1, n = [1, 0])
            get_neighbor_face(ii, jj, 3, u_plus);
            for (int k = 0; k < 3; k++) {
                double uM = cell->face_right[k];
                double uP = u_plus[k];
                double fn_int = flux_f(uM);
                double fn_num = lf_flux(uM,uP,1,0);
                double diff = fn_int - fn_num;
                for (int i = 0; i < N9; i++) 
                    surf_res[i] += (dy / 2.0) * weights1D[k] * diff * face_interp_R[k][i];
            }

            /* --- 步骤 C: 乘以质量矩阵之逆并合并 --- */
            // 直接加上边界残差向量，杜绝循环内的标量化
            for (int i = 0; i < N9; i++) {
                RHS[idx][i] += surf_res[i];
            }

            double RHS_temp[N9] = {0};
            for (int i = 0; i < N9; i++) {
                for(int col = 0; col < N9; col++) {
                    RHS_temp[i] += (M_inv[i][col] * RHS[idx][col]);
                }
            }
            
            // 必须除以参考雅可比行列式 J
            for (int i = 0; i < N9; i++) {
                RHS[idx][i] = RHS_temp[i] / J;
            }
        }
    }
}


/* ============================================================
 * 时间步长（CFL 条件）
 * ============================================================ */
static double compute_dt(void) {
    double speed = max_speed();
    double h     = (dx < dy) ? dx : dy;
    return CFL * h / speed;
}

/* ============================================================
 * 三阶 SSP-RK3 时间推进
 * ============================================================ */
static double U0[Nx * Ny][N9]; /* U^n  暂存 */
static double U1[Nx * Ny][N9]; /* U^(1) */
static double U2[Nx * Ny][N9]; /* U^(2) */
static double RHS_buf[Nx * Ny][N9];

void rk3_step(double dt,int nit) {
    int total = Nx * Ny;

    /* 保存 U^n */
    for (int i = 0; i < total; i++)
        memcpy(U0[i], Mesh[i].U, N9 * sizeof(double));

    /* Stage 1: U^(1) = U^n + dt * L(U^n) */
    compute_rhs(RHS_buf);
   

    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            U1[i][q] = U0[i][q] + dt * RHS_buf[i][q];
    for (int i = 0; i < total; i++)
        memcpy(Mesh[i].U, U1[i], N9 * sizeof(double));
apply_limiter();
    /* Stage 2: U^(2) = 3/4 U^n + 1/4 (U^(1) + dt * L(U^(1))) */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            U2[i][q] = 0.75 * U0[i][q]
                     + 0.25 * (U1[i][q] + dt * RHS_buf[i][q]);
    for (int i = 0; i < total; i++)
        memcpy(Mesh[i].U, U2[i], N9 * sizeof(double));
apply_limiter();
    /* Stage 3: U^{n+1} = 1/3 U^n + 2/3 (U^(2) + dt * L(U^(2))) */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            Mesh[i].U[q] = (1.0 / 3.0) * U0[i][q]
                         + (2.0 / 3.0) * (U2[i][q] + dt * RHS_buf[i][q]);
apply_limiter();

}


void output_results(double t) {
    char filename[50];
    sprintf(filename, "result.dat");
    FILE *fp = fopen(filename, "w");
    
    if (fp == NULL) {
        printf("Error: Could not open file for writing!\n");
        return;
    }
    fprintf(fp, "VARIABLES = \"X\", \"Y\", \"U\"\n");

    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
       
            int idx = jj * Nx + ii;
            Element *cell = &Mesh[idx];

            // 计算当前单元的中心坐标
            double xc = (ii + 0.5) * dx;
            double yc = (jj + 0.5) * dy;

            // 遍历单元内的 9 个节点 (3x3 布局)
            for (int i = 0; i < 3; i++) {       // y 方向索引
                for (int j = 0; j < 3; j++) {   // x 方向索引
                    int node_idx = i * 3 + j;
                    
                    // 将参考坐标 (-1, 0, 1) 映射到物理坐标
                    // 映射公式: x = xc + (dx/2) * xi
                    double x_phys = xc + (dx / 2.0) * nodes1D[j];
                    double y_phys = yc + (dy / 2.0) * nodes1D[i];
            
                    fprintf(fp, "%lf %lf %lf \n", x_phys, y_phys, cell->U[node_idx]);
                }
            }
        }
    }

    fclose(fp);
    printf("Saved results to %s\n", filename);
}
/* ============================================================
 * 主函数
 * ============================================================ */
int main(void) {
    init_quadrature();
    precompute_matrices();
    init_condition();

    double t  = 0.0;
    int   nit = 0;

    printf("%-10s  %-14s  %-14s\n", "Step", "Time", "L_error");
    printf("--------------------------------------------\n");

    while (t < T_END) {
        double dt = compute_dt();
        if (t + dt > T_END) dt = T_END - t;

        rk3_step(dt,nit);
        t  += dt;
        nit++;

        if (nit % 50 == 0 || t >= T_END) {
     
            printf("%-10d  %-14.6e  \n", nit, t);
        }
    }
output_results(T_END);
    return 0;
}