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
#define T_END 1
#define deltat 0.3
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
static const double ADV_A = 1.0; /* x 方向波速 */
static const double ADV_B = 1.0; /* y 方向波速 */

/* x 方向通量 f(u) */
static inline double flux_f(double u) {
    return ADV_A * u;
}

/* y 方向通量 g(u) */
static inline double flux_g(double u) {
    return ADV_B * u;
}

/* ============================================================
 * 最大波速（用于 LF 通量稳定性参数 alpha）
 * ============================================================ */
static inline double max_speed_x(void) { return fabs(ADV_A); }
static inline double max_speed_y(void) { return fabs(ADV_B); }

/* ============================================================
 * LF（Lax-Friedrichs）数值通量
 * ============================================================ */
static inline double lf_flux_x(double u_minus, double u_plus) {
    double alpha = max_speed_x();
    return 0.5 * (flux_f(u_minus) + flux_f(u_plus))
           - 0.5 * alpha * (u_plus - u_minus);
}

static inline double lf_flux_y(double u_minus, double u_plus) {
    double alpha = max_speed_y();
    return 0.5 * (flux_g(u_minus) + flux_g(u_plus))
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
            /* 单元左下角物理坐标 */
            double x0 = xL + ii * dx;
            double y0 = yL + jj * dy;
            for (int q = 0; q < 9; q++) {
                /* 参考坐标 -> 物理坐标 */
                double x = x0 + 0.5 * dx * (r_quad[q] + 1.0);
                double y = y0 + 0.5 * dy * (s_quad[q] + 1.0);
                cell->U[q] = sin(2.0 * M_PI * (x-deltat)) * sin(2.0 * M_PI * (y-deltat));
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
 * Ghost Cell 边界条件（周期性）

 * 若需 Dirichlet 边界，将对应 ghost face 值替换为边界值即可。
 * ============================================================ */
void apply_ghost_cells(void) {
    /* 底/顶 Ghost */
    for (int i = 0; i < Nx; i++) {
        for (int k = 0; k < 3; k++) {
            Ghost_bottom[i].face_top[k]    = Mesh[(Ny - 1) * Nx + i].face_top[k];
            Ghost_top[i].face_bottom[k]    = Mesh[0 * Nx + i].face_bottom[k];
        }
    }
    /* 左/右 Ghost */
    for (int j = 0; j < Ny; j++) {
        for (int k = 0; k < 3; k++) {
            Ghost_left[j].face_right[k]    = Mesh[j * Nx + (Nx - 1)].face_right[k];
            Ghost_right[j].face_left[k]    = Mesh[j * Nx + 0].face_left[k];
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
                double fn_num = 0.5 * (-flux_g(uM) - flux_g(uP)) - 0.5 * max_speed_y() * (uP - uM);
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
                double fn_num = 0.5 * (flux_g(uM) + flux_g(uP)) - 0.5 * max_speed_y() * (uP - uM);
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
                double fn_num = 0.5 * (-flux_f(uM) - flux_f(uP)) - 0.5 * max_speed_x() * (uP - uM);
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
                double fn_num = 0.5 * (flux_f(uM) + flux_f(uP)) - 0.5 * max_speed_x() * (uP - uM);
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
    double speed = max_speed_x() + max_speed_y();
    double h     = (dx < dy) ? dx : dy;
    return CFL * h / speed;
}

/* ============================================================
 * 三阶 SSP-RK3 时间推进
 *
 *   U^(1) = U^n + dt * L(U^n)
 *   U^(2) = (3/4) U^n + (1/4)[U^(1) + dt * L(U^(1))]
 *   U^n+1 = (1/3) U^n + (2/3)[U^(2) + dt * L(U^(2))]
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
/*
    for (int i = 0; i < Nx * Ny; i++) {
        // 打印单元索引
        printf("Cell %-4d: ", i);
        for (int j = 0; j < N9; j++) {
            // 使用 %.4e (科学计数法) 既精确又整齐
            printf("%12.4e ", RHS_buf[i][j]);
        }
        putchar('\n');
    }
*/
   

    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            U1[i][q] = U0[i][q] + dt * RHS_buf[i][q];
    for (int i = 0; i < total; i++)
        memcpy(Mesh[i].U, U1[i], N9 * sizeof(double));

    /* Stage 2: U^(2) = 3/4 U^n + 1/4 (U^(1) + dt * L(U^(1))) */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            U2[i][q] = 0.75 * U0[i][q]
                     + 0.25 * (U1[i][q] + dt * RHS_buf[i][q]);
    for (int i = 0; i < total; i++)
        memcpy(Mesh[i].U, U2[i], N9 * sizeof(double));

    /* Stage 3: U^{n+1} = 1/3 U^n + 2/3 (U^(2) + dt * L(U^(2))) */
    compute_rhs(RHS_buf);
    for (int i = 0; i < total; i++)
        for (int q = 0; q < N9; q++)
            Mesh[i].U[q] = (1.0 / 3.0) * U0[i][q]
                         + (2.0 / 3.0) * (U2[i][q] + dt * RHS_buf[i][q]);


}

/* ============================================================
 * L2 误差（精确解已知时用于验证）
 * ============================================================ */
double compute_l2_error(double t) {
    double err2 = 0.0;
    for (int jj = 0; jj < Ny; jj++) {
        for (int ii = 0; ii < Nx; ii++) {
            Element *cell = &Mesh[jj * Nx + ii];
            double x0 = xL + ii * dx;
            double y0 = yL + jj * dy;
            /* 2D 张量积 Gauss 积分 */
            for (int q = 0; q < 9; q++) {
                double x    = x0 + 0.5 * dx * (r_quad[q] + 1.0);
                double y    = y0 + 0.5 * dy * (s_quad[q] + 1.0);
                /* 精确解（周期边界）：u_ex = sin(2pi(x-at)) * sin(2pi(y-bt)) */
                double uex  = sin(2.0 * M_PI * (x - ADV_A * (t+deltat)))
                            * sin(2.0 * M_PI * (y - ADV_B * (t+deltat)));
                double diff = cell->U[q] - uex;
                /* Jacobian = (dx/2)*(dy/2) */
               err2 += w_quad[q] * diff * diff * (dx / 2.0) * (dy / 2.0);
             
            }
        }
    }
     err2=sqrt(err2);
    return (err2);
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
                       double uex  = sin(2.0 * M_PI * (x_phys - ADV_A * (t+deltat)))
                            * sin(2.0 * M_PI * (y_phys - ADV_B * (t+deltat)));
                    fprintf(fp, "%lf %lf %lf %lf\n", x_phys, y_phys, cell->U[node_idx],uex);
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
            double err = compute_l2_error(t);
            printf("%-10d  %-14.6e  %-14.6e\n", nit, t, err);
        }
    }
output_results(T_END);
    return 0;
}