% 读取 sedov_1d.dat 文件并绘制结果
% 运行前请确保 sedov_1d.dat 与此脚本位于同一目录下

% 使用 readmatrix 读取数据，跳过第一行表头
filename = 'sedov_1d.dat';
if ~isfile(filename)
    error('未找到文件 %s，请检查 C 程序是否成功运行并生成了该文件。', filename);
end

data = readmatrix(filename, 'NumHeaderLines', 1);

% 提取各个物理量
x   = data(:, 1);  % 空间坐标 x
rho = data(:, 2);  % 密度
u   = data(:, 3);  % 速度
p   = data(:, 4);  % 压力

% 创建绘图窗口
figure;


plot(x, log10(rho), '-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'b');
xlabel('x');
ylabel('\rho');
title('Density');
grid on;


