%% DG 结果可视化 - 坐标对齐重构（黑白）
clear; clc;

% 1. 加载数据
data = readmatrix('320_160.dat', 'FileType', 'text', 'NumHeaderLines', 1);
x_raw = data(:, 1); 
y_raw = data(:, 2); 
u_raw = data(:, 3); 



% 1. 提取第 6 列数据
column6 = data(:, 6)+1e-13;

% 2. 定义筛选条件
% condition1: 负数或零 (实数范围内)
% condition2: 包含虚部的复数
% condition3: 异常值 NaN
is_negative_or_zero = (column6 <= 0);
is_complex = ~isreal(column6);
is_not_a_number = isnan(column6);

% 3. 合并所有“非正实数”或“异常”的逻辑索引
% 使用逻辑“或”操作符 |
target_mask = is_negative_or_zero | is_complex | is_not_a_number;

% 4. 提取并打印结果
invalid_data = column6(target_mask);

disp('找到的负数、0 或非实数如下：');
disp(invalid_data);




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
rho_log=log10(U)/log10(2.5);
% 4. 检查是否有空洞（若 DG 节点分布不均，矩阵中会有 NaN）
if any(isnan(U(:)))
    warning('网格中存在未定义的点，等值线可能会断裂。');
end


[X, Y] = meshgrid(ux, uy);



figure('Color', 'w');

% 使用填充等值线图 (contourf) 通常比普通 contour 在对数尺度下更好看
levels = linspace(min(rho_log(:)), max(rho_log(:)), 1000);
contourf(X, Y, rho_log, levels, 'LineColor', 'none');

% 1. 定义纯彩虹色带的 5 个基础颜色节点：蓝 -> 青 -> 绿 -> 黄 -> 红
base_colors = [
    0, 0, 0.6;  % Deep Blue (最低)
    0, 0.5, 1;   % Light Blue (过渡)
    0, 1, 1;     % Cyan (过渡)
    1, 1, 0;     % Yellow (中间 - 节点 1)
    1, 1, 0;     % Yellow (中间 - 节点 2，重复以扩大区域)
    1, 0, 0      % Red (最高)
];

% 2. 定义这 5 个节点在 0 到 1 区间上的相对位置
nodes = linspace(0, 1, 6);

% 3. 生成具有 256 个层级的平滑自定义色带
query_points = linspace(0, 1, 256)';
custom_colormap = interp1(nodes, base_colors, query_points);

% 4. 应用自定义色带，并严格锁死数据映射范围
colormap(custom_colormap);
clim([-2, 3]); % 如果 MATLAB 版本较老，请使用 caxis([-2, 3]);

% 5. 设置颜色条刻度
cb = colorbar;
cb.Ticks = -2:0.5:3;

axis equal;
xlabel('x'); ylabel('y');
shading interp; % 平滑显示

disp(max(max(U)));

