%% DG 结果可视化脚本
clear; clc;

% 1. 加载数据
filename = 'result.dat'; % 对应你步数 300 的文件
data = readmatrix(filename, 'FileType', 'text', 'NumHeaderLines', 1);

x = data(:, 1); % 第一列: X 坐标
y = data(:, 2); % 第二列: Y 坐标
u = data(:, 3); % 第三列: 变量 U

% 2. 绘制散点云图 (初步检查)
figure(1);
scatter3(x, y, u, 20, u, 'filled');
colorbar;
xlabel('X'); ylabel('Y'); zlabel('U');
title('DG 离散节点分布图');
view(2); % 切换到俯视图

% 3. 绘制平滑曲面图 (精美可视化)
% 创建插值网格
[xq, yq] = meshgrid(min(x):0.01:max(x), min(y):0.01:max(y));
F = scatteredInterpolant(x, y, u, 'linear');
uq = F(xq, yq);

figure(2);
surf(xq, yq, uq, 'EdgeColor', 'none');
shading interp; % 平滑着色
colormap jet;
colorbar;
xlabel('X'); ylabel('Y'); zlabel('U');
title('2D DG 数值解可视化');

camlight;

