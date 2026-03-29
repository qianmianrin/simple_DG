#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 20             // X方向网格数
#define NY 20             // Y方向网格数
#define LX 2.0            // X方向区域长度
#define LY 2.0            // Y方向区域长度
#define N_VAR 4           // 变量数 (rho, rhou, rhov, E)
#define N_DEG 2           // 多项式最高阶数
#define N_BASIS 9         // Q2空间基函数数量 (3x3)
#define N_GP 3            // 一维高斯点数
#define GAMMA 1.4         // 气体常数
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// 3点高斯-勒让德求积节点与权重
const double gp[3] = {-0.774596669241483, 0.0, 0.774596669241483};
const double gw[3] = {0.555555555555555, 0.888888888888888, 0.555555555555555};

// 单元结构体
typedef struct {
    double U[N_VAR][N_BASIS];
    double RHS[N_VAR][N_BASIS];
} Element;

Element mesh[NX][NY];
Element mesh_U1[NX][NY];
Element mesh_U2[NX][NY];
double dx, dy;
double mass_inv[N_BASIS];

// 勒让德多项式及其导数
double L(int m, double x) {
    if (m == 0) return 1.0;
    if (m == 1) return x;
    return 0.5 * (3.0 * x * x - 1.0);
}

double dL(int m, double x) {
    if (m == 0) return 0.0;
    if (m == 1) return 1.0;
    return 3.0 * x;
}

// 获取某一状态的点值
void get_state(double U_coeffs[N_VAR][N_BASIS], double xi, double eta, double U[N_VAR]) {
    for (int v = 0; v < N_VAR; v++) U[v] = 0.0;
    int k = 0;
    for (int my = 0; my <= N_DEG; my++) {
        for (int mx = 0; mx <= N_DEG; mx++) {
            double phi = L(mx, xi) * L(my, eta);
            for (int v = 0; v < N_VAR; v++) U[v] += U_coeffs[v][k] * phi;
            k++;
        }
    }
}

// 物理通量计算
void get_flux(double U[N_VAR], double F[N_VAR], double G[N_VAR]) {
    double rho = U[0], u = U[1]/rho, v = U[2]/rho, E = U[3];
    double p = (GAMMA - 1.0) * (E - 0.5 * rho * (u*u + v*v));
    F[0] = rho * u;          G[0] = rho * v;
    F[1] = rho * u * u + p;  G[1] = rho * u * v;
    F[2] = rho * u * v;      G[2] = rho * v * v + p;
    F[3] = (E + p) * u;      G[3] = (E + p) * v;
}

// 声速及最大特征值
double max_eig(double U[N_VAR], int dir) {
    double rho = U[0], u = U[1]/rho, v = U[2]/rho, E = U[3];
    double p = (GAMMA - 1.0) * (E - 0.5 * rho * (u*u + v*v));
    double c = sqrt(GAMMA * p / rho);
    return (dir == 0) ? (fabs(u) + c) : (fabs(v) + c);
}

// Local Lax-Friedrichs 数值通量
void llf_flux(double UL[N_VAR], double UR[N_VAR], int dir, double flux[N_VAR]) {
    double FL[N_VAR], GL[N_VAR], FR[N_VAR], GR[N_VAR];
    get_flux(UL, FL, GL);
    get_flux(UR, FR, GR);
    double lambda = fmax(max_eig(UL, dir), max_eig(UR, dir));
    
    for (int v = 0; v < N_VAR; v++) {
        if (dir == 0) flux[v] = 0.5 * (FL[v] + FR[v] - lambda * (UR[v] - UL[v]));
        else          flux[v] = 0.5 * (GL[v] + GR[v] - lambda * (UR[v] - UL[v]));
    }
}


void init_cond() {
    dx = LX / NX; dy = LY / NY;
    // 预计算质量矩阵逆
    int k = 0;
    for (int my = 0; my <= N_DEG; my++) {
        for (int mx = 0; mx <= N_DEG; mx++) {
            double mass = (dx * dy / 4.0) * (2.0 / (2*mx + 1)) * (2.0 / (2*my + 1));
            mass_inv[k++] = 1.0 / mass;
        }
    }

    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            double cx = i * dx + dx / 2.0;
            double cy = j * dy + dy / 2.0;
            
            for(int v=0; v<N_VAR; v++) 
                for(int b=0; b<N_BASIS; b++) mesh[i][j].U[v][b] = 0.0;

            for (int gy = 0; gy < N_GP; gy++) {
                for (int gx = 0; gx < N_GP; gx++) {
                    double xi = gp[gx], yi = gp[gy];
                    double weight = gw[gx] * gw[gy] * (dx * dy / 4.0);
                    double x = cx + xi * dx / 2.0;
                    double y = cy + yi * dy / 2.0;
                    
                    //初始条件
                    double rho = 1.0 + 0.5 * sin(M_PI * (x + y));
                    double u = 1.0, v_vel = 1.0, p = 1.0;
                    double E = p / (GAMMA - 1.0) + 0.5 * rho * (u*u + v_vel*v_vel);
                    double U_exact[N_VAR] = {rho, rho*u, rho*v_vel, E};

                    int k_b = 0;
                    for (int my = 0; my <= N_DEG; my++) {
                        for (int mx = 0; mx <= N_DEG; mx++) {
                            double phi = L(mx, xi) * L(my, yi);
                            for(int v=0; v<N_VAR; v++) mesh[i][j].U[v][k_b] += U_exact[v] * phi * weight;
                            k_b++;
                        }
                    }
                }
            }
            // 乘以质量矩阵逆
            for(int v=0; v<N_VAR; v++) 
                for(int b=0; b<N_BASIS; b++) mesh[i][j].U[v][b] *= mass_inv[b];
        }
    }
}

// 计算当前状态的空间残差
void compute_rhs(Element m[NX][NY]) {
    // 1. 初始化 RHS 为 0
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            for (int v = 0; v < N_VAR; v++)
                for (int b = 0; b < N_BASIS; b++)
                    m[i][j].RHS[v][b] = 0.0;

    // 2. 体积分 (Volume Integral)
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int gy = 0; gy < N_GP; gy++) {
                for (int gx = 0; gx < N_GP; gx++) {
                    double xi = gp[gx], eta = gp[gy];
                    double w = gw[gx] * gw[gy] * (dx * dy / 4.0);
                    
                    double U_val[N_VAR], F[N_VAR], G[N_VAR];
                    get_state(m[i][j].U, xi, eta, U_val);
                    get_flux(U_val, F, G);

                    int k = 0;
                    for (int my = 0; my <= N_DEG; my++) {
                        for (int mx = 0; mx <= N_DEG; mx++) {
                            double dphi_dx = dL(mx, xi) * L(my, eta) * (2.0 / dx);
                            double dphi_dy = L(mx, xi) * dL(my, eta) * (2.0 / dy);
                            for (int v = 0; v < N_VAR; v++) {
                                m[i][j].RHS[v][k] += w * (F[v] * dphi_dx + G[v] * dphi_dy);
                            }
                            k++;
                        }
                    }
                }
            }
        }
    }

    // 3. 面积分 (Surface Integral) - X 方向边界 (利用周期性边界)
    for (int i = 0; i < NX; i++) {
        int i_right = (i + 1) % NX;
        for (int j = 0; j < NY; j++) {
            for (int g = 0; g < N_GP; g++) {
                double eta = gp[g];
                double w = gw[g] * (dy / 2.0);
                
                double UL[N_VAR], UR[N_VAR], flux[N_VAR];
                get_state(m[i][j].U, 1.0, eta, UL);          // 左单元右边界
                get_state(m[i_right][j].U, -1.0, eta, UR);   // 右单元左边界
                llf_flux(UL, UR, 0, flux);

                int k = 0;
                for (int my = 0; my <= N_DEG; my++) {
                    for (int mx = 0; mx <= N_DEG; mx++) {
                        double phi_L = L(mx, 1.0) * L(my, eta);
                        double phi_R = L(mx, -1.0) * L(my, eta);
                        for (int v = 0; v < N_VAR; v++) {
                            m[i][j].RHS[v][k]       -= w * flux[v] * phi_L;
                            m[i_right][j].RHS[v][k] += w * flux[v] * phi_R;
                        }
                        k++;
                    }
                }
            }
        }
    }

    // 4. 面积分 (Surface Integral) - Y 方向边界
    for (int j = 0; j < NY; j++) {
        int j_top = (j + 1) % NY;
        for (int i = 0; i < NX; i++) {
            for (int g = 0; g < N_GP; g++) {
                double xi = gp[g];
                double w = gw[g] * (dx / 2.0);
                
                double UB[N_VAR], UT[N_VAR], flux[N_VAR];
                get_state(m[i][j].U, xi, 1.0, UB);           // 下单元上边界
                get_state(m[i][j_top].U, xi, -1.0, UT);      // 上单元下边界
                llf_flux(UB, UT, 1, flux);

                int k = 0;
                for (int my = 0; my <= N_DEG; my++) {
                    for (int mx = 0; mx <= N_DEG; mx++) {
                        double phi_B = L(mx, xi) * L(my, 1.0);
                        double phi_T = L(mx, xi) * L(my, -1.0);
                        for (int v = 0; v < N_VAR; v++) {
                            m[i][j].RHS[v][k]       -= w * flux[v] * phi_B;
                            m[i][j_top].RHS[v][k]   += w * flux[v] * phi_T;
                        }
                        k++;
                    }
                }
            }
        }
    }
}

// 输出密度场的单元平均值到文件
void output_results(int nx, int ny) {
    char grid_fn[64], dens_fn[64];
    sprintf(grid_fn, "grid_Nx%d_Ny%d.txt", nx, ny);
    sprintf(dens_fn, "density_Nx%d_Ny%d.txt", nx, ny);

    // 1. 输出网格文件 (xc, yc)
    FILE *fg = fopen(grid_fn, "w");
    fprintf(fg, "%d %d\n", nx, ny);
    for (int i = 0; i < nx; i++) fprintf(fg, "%e ", (i + 0.5) * dx);
    fprintf(fg, "\n");
    for (int j = 0; j < ny; j++) fprintf(fg, "%e ", (j + 0.5) * dy);
    fclose(fg);

    // 2. 输出密度文件 (尺寸 + 矩阵)
    FILE *fd = fopen(dens_fn, "w");
    fprintf(fd, "%d %d\n", nx, ny);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fprintf(fd, "%e ", mesh[i][j].U[0][0]); // U[0][0] 是 Q2 的单元平均项
        }
        fprintf(fd, "\n");
    }
    fclose(fd);
    printf("Results saved to %s and %s\n", grid_fn, dens_fn);
}
int main() {
    init_cond();
    
    double t = 0.0, t_end = 2.0;
    double cfl = 0.1;
    // 根据声速估算稳定时间步长 (Max wave speed approx 2.2 for this IC)
    double dt = cfl * dx / (2.2 * (2 * N_DEG + 1)); 

    printf("Starting DG Q2 simulation. dt = %f\n", dt);
    
    int step = 0;
    while (t < t_end) {
        if (t + dt > t_end) dt = t_end - t;

        // --- SSP-RK3 Stage 1 ---
        compute_rhs(mesh);
        for (int i=0; i<NX; i++) for (int j=0; j<NY; j++) 
            for (int v=0; v<N_VAR; v++) for (int b=0; b<N_BASIS; b++)
                mesh_U1[i][j].U[v][b] = mesh[i][j].U[v][b] + dt * mesh[i][j].RHS[v][b] * mass_inv[b];

        // --- SSP-RK3 Stage 2 ---
        compute_rhs(mesh_U1);
        for (int i=0; i<NX; i++) for (int j=0; j<NY; j++) 
            for (int v=0; v<N_VAR; v++) for (int b=0; b<N_BASIS; b++)
                mesh_U2[i][j].U[v][b] = 0.75 * mesh[i][j].U[v][b] + 0.25 * mesh_U1[i][j].U[v][b] + 
                                        0.25 * dt * mesh_U1[i][j].RHS[v][b] * mass_inv[b];

        // --- SSP-RK3 Stage 3 ---
        compute_rhs(mesh_U2);
        for (int i=0; i<NX; i++) for (int j=0; j<NY; j++) 
            for (int v=0; v<N_VAR; v++) for (int b=0; b<N_BASIS; b++)
                mesh[i][j].U[v][b] = 1.0/3.0 * mesh[i][j].U[v][b] + 2.0/3.0 * mesh_U2[i][j].U[v][b] + 
                                     2.0/3.0 * dt * mesh_U2[i][j].RHS[v][b] * mass_inv[b];

        t += dt;
        step++;
        if (step % 100 == 0) printf("Step %d, Time %f\n", step, t);
    }

    output_results(NX,NY);
    return 0;
}