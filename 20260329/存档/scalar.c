#include<stdio.h>
#include <stdlib.h>
#include <math.h>
//基本全局参数
#define Nx 30
#define Ny 30
#define xL 0
#define xR 1
#define yL 0
#define yR 1
#define dx (xR-xL)/Nx
#define dy (yR-yL)/Ny

//网格
typedef struct {
    double U[9];
    double face_top[3];
    double face_bottom[3];
    double face_left[3];
    double face_right[3];
} Element;
//三阶Gauss-Legend积分点，共九个

Element Mesh[Nx*Ny];//第i，j个单元对应Mesh[j*Nx+i]

static const double nodes1D[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
static const double weights1D[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

// 二维张量积节点 (r, s) 共 9 个
double r_quad[9], s_quad[9], w_quad[9];

void init_quadrature() {
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            int idx = j * 3 + i;
            r_quad[idx] = nodes1D[i];
            s_quad[idx] = nodes1D[j];
            w_quad[idx] = weights1D[i] * weights1D[j];
        }
    }
}
//基函数，导数

static inline double lagrange_1d(int i, double x) {
    double val = 1.0;
    for (int j = 0; j < 3; j++) {
        if (i == j) continue;
        val *= (x - nodes1D[j]) / (nodes1D[i] - nodes1D[j]);
    }
    return val;
}

// 1D Lagrange 导数 dL_i/dx
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




// 预计算矩阵全局变量 (9x9 对应 2D, 3x3 对应 1D)
double M_inv[9][9];      // 逆质量矩阵
double Dr[9][9];         // r 方向微分矩阵 (M^-1 * Sr)
double Ds[9][9];         // s 方向微分矩阵 (M^-1 * Ss)

// 边界插值矩阵：从 9 个内部点插值到 4 条边的 3 个 Gauss 点
double face_interp_T[3][9];
double face_interp_B[3][9];
double face_interp_L[3][9];
double face_interp_R[3][9];

void precompute_matrices() {
    double M[9][9] = {0};
    double Sr[9][9] = {0}; // Stiffness matrix in r
    double Ss[9][9] = {0}; // Stiffness matrix in s


    // 组装体积分矩阵 (Mass and Stiffness)
    // 索引映射说明：idx = m*3 + n 对应 (node_n, node_m)
    for (int q = 0; q < 9; q++) { // 遍历 9 个积分点
        double r = r_quad[q];
        double s = s_quad[q];
        double w = w_quad[q];

        for (int i = 0; i < 9; i++) {
            int ni = i % 3; int mi = i / 3;
            double phi_i = lagrange_1d(ni, r) * lagrange_1d(mi, s);
            
            for (int j = 0; j < 9; j++) {
                int nj = j % 3; int mj = j / 3;
                double phi_j = lagrange_1d(nj, r) * lagrange_1d(mj, s);
                double dphi_j_dr = d_lagrange_1d(nj, r) * lagrange_1d(mj, s);
                double dphi_j_ds = lagrange_1d(nj, r) * d_lagrange_1d(mj, s);

                M[i][j]  += w * phi_i * phi_j;
                Sr[i][j] += w * phi_i * dphi_j_dr;
                Ss[i][j] += w * phi_i * dphi_j_ds;
            }
        }
    }

    // 计算面插值矩阵 
    for (int k = 0; k < 3; k++) { // 边界上的 3 个 Gauss 点
        double node_k = nodes1D[k];
        for (int j = 0; j < 9; j++) {
            int n = j % 3; int m = j / 3;
            face_interp_T[k][j] = lagrange_1d(n, node_k) * lagrange_1d(m,  1.0);
            face_interp_B[k][j] = lagrange_1d(n, node_k) * lagrange_1d(m, -1.0);
            face_interp_L[k][j] = lagrange_1d(n, -1.0)   * lagrange_1d(m, node_k);
            face_interp_R[k][j] = lagrange_1d(n,  1.0)   * lagrange_1d(m, node_k);
        }
    }

    // 4. 矩阵求逆及微分矩阵计算 (此处简记，实际应用中建议使用 LU 分解)
    // 对于 9x9 矩阵，可以使用简单的辅助函数或硬编码求逆
    // 得到 M_inv 后：
    // Dr = M_inv * Sr;
    // Ds = M_inv * Ss;
}

//到边界点的计算，使用预计算矩阵进行边界插值
void boundary_value() {
    for (int idx = 0; idx < Nx * Ny; idx++) {
        Element *cell = &Mesh[idx];

        // 遍历边界上的 3 个 1D Gauss 点 (index k)
        for (int k = 0; k < 3; k++) {
            double val_top = 0.0;
            double val_bottom = 0.0;
            double val_left = 0.0;
            double val_right = 0.0;

            // 核心优化：利用预计算的插值矩阵进行线性组合
            // j 遍历单元内的 9 个 2D 节点
            for (int j = 0; j < 9; j++) {
                double u_val = cell->U[j];
                
                val_top    += face_interp_T[k][j] * u_val;
                val_bottom += face_interp_B[k][j] * u_val;
                val_left   += face_interp_L[k][j] * u_val;
                val_right  += face_interp_R[k][j] * u_val;
            }

            // 存储结果
            cell->face_top[k]    = val_top;
            cell->face_bottom[k] = val_bottom;
            cell->face_left[k]   = val_left;
            cell->face_right[k]  = val_right;
        }
    }
}

void main(){
    init_quadrature();
}