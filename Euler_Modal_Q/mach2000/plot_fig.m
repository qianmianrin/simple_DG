%% DG 结果可视化 - 坐标对齐重构（黑白）
clear; clc;

% 1. 加载数据
data = readmatrix('result.dat', 'FileType', 'text', 'NumHeaderLines', 1);
x_raw = data(:, 1); 
y_raw = data(:, 2); 
u_raw = data(:, 3); 

% 2. 确定全局坐标范围和分辨率
% 我们依然使用 unique 获取坐标轴线，但仅用于建立“刻度”
ux = unique(x_raw); 
uy = unique(y_raw);
nx = length(ux); 
ny = length(uy);

% 3. 建立坐标到索引的映射 (核心：解决序关系破坏问题)
% 预分配一个空矩阵
U = NaN(ny, nx); 

% 遍历原始数据，将每个 u 放到它物理坐标对应的矩阵位置上
% ismemberloc 是为了找到原始坐标在 ux, uy 向量中的位置索引
[~, x_idx] = ismember(x_raw, ux);
[~, y_idx] = ismember(y_raw, uy);

% 线性索引转换：将 (y_idx, x_idx) 映射到矩阵 U
% sub2ind 保证了物理上的 (x,y) 正确对应到矩阵的行和列
linear_idx = sub2ind([ny, nx], y_idx, x_idx);
U(linear_idx) = u_raw;
rho_log=log10(U);
% 4. 检查是否有空洞（若 DG 节点分布不均，矩阵中会有 NaN）
if any(isnan(U(:)))
    warning('网格中存在未定义的点，等值线可能会断裂。');
end


[X, Y] = meshgrid(ux, uy);



figure('Color', 'w');

% 使用填充等值线图 (contourf) 通常比普通 contour 在对数尺度下更好看
levels = linspace(min(rho_log(:)), max(rho_log(:)), 100);
contourf(X, Y, rho_log, levels, 'LineColor', 'none');

% 设置配色方案
colormap(jet(256));
colorbar;

axis equal;
xlabel('x'); ylabel('y');
shading interp; % 平滑显示
