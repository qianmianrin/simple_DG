function plot_density_euler(Nx, Ny)
% plot_density_euler - 针对二维 Euler DG 程序定制的绘图脚本
% 使用方法: plot_density_euler(20, 20)

if nargin < 2
    Nx = 20;
    Ny = 20;
end

%% 1. 读取网格中心点坐标
grid_file = sprintf('grid_Nx%d_Ny%d.txt', Nx, Ny);
if ~exist(grid_file, 'file'), error('找不到网格文件: %s', grid_file); end

fid = fopen(grid_file, 'r');
grid_dims = fscanf(fid, '%d %d', 2);
xc = fscanf(fid, '%e', grid_dims(1));
yc = fscanf(fid, '%e', grid_dims(2));
fclose(fid);

%% 2. 读取密度场数据
density_file = sprintf('density_Nx%d_Ny%d.txt', Nx, Ny);
if ~exist(density_file, 'file'), error('找不到密度文件: %s', density_file); end

fid = fopen(density_file, 'r');
data_dims = fscanf(fid, '%d %d', 2);
% 注意：fscanf 读取时是按列填充，由于 C 写入是行优先，
% 我们读入 [Nx, Ny] 的矩阵后，rho_data(i,j) 对应 C 里的 mesh[i][j]
rho_data = fscanf(fid, '%e', [data_dims(1), data_dims(2)]);
fclose(fid);

% 为了配合 meshgrid(xc, yc) 的坐标映射，我们需要转置，使得维度匹配 [length(yc), length(xc)]
rho = rho_data'; 

%% 3. 绘图与可视化
[X, Y] = meshgrid(xc, yc);

figure;
% 使用 100 层等值线，并开启填充
contour(X, Y, rho, linspace(min(min(rho)),max(max(rho)),50), 'k','LineWidth', 0.8);


axis equal;
xlabel('Position x'); ylabel('Position y');
title(sprintf('Density Field (Nx=%d, Ny=%d) - Smooth Solver', Nx, Ny));

%% 4. 显示统计信息
fprintf('========================================\n');
fprintf('Euler DG 2D 统计 (Q2 Polynomials):\n');
fprintf('----------------------------------------\n');
fprintf('  数据规模: %d x %d\n', data_dims(1), data_dims(2));
fprintf('  密度极值: [%.6e, %.6e]\n', min(rho(:)), max(rho(:)));
fprintf('  总质量近似: %.6f\n', sum(rho(:)) * (xc(2)-xc(1)) * (yc(2)-yc(1)));
fprintf('========================================\n');
end