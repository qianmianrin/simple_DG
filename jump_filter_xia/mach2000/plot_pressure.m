% 1. 建立坐标映射（假设 x_raw, y_raw, column6 已经准备好）
U = NaN(ny, nx); 
[~, x_idx] = ismember(x_raw, ux);
[~, y_idx] = ismember(y_raw, uy);
linear_idx = sub2ind([ny, nx], y_idx, x_idx);
U(linear_idx) = column6; % 填入压强数据，不要加 1e-13

% 2. 剔除负数和0，防止复数
U_valid_mask = (U > 0);
U_plot = U; 
U_plot(~U_valid_mask) = NaN; 
rho_log = log10(U_plot);

% 3. 绘图
[X, Y] = meshgrid(ux, uy);
figure('Color', 'w');

% 直接让 MATLAB 画 100 层等值线
contourf(X, Y, rho_log, 1000, 'LineColor', 'none'); 

% 使用 jet 色系
colormap('jet'); 

% ================== 核心修正 ==================
% 物理背景压强是 0.4127，取对数约为 -0.38。
% 我们把颜色条的下限强制锁死在 -0.5，这样 -0.38 就会落在深蓝区域。
% 那些小于 -0.5 的微小数值误差依然存在，但它们全都会被统一涂成深蓝色，不会再污染整个流场的配色。
physical_min_log = -0.5; 
actual_max_log = max(rho_log(U_valid_mask));

% 应用颜色范围限制 (如果是 R2022a 之前的 MATLAB，请把 clim 换成 caxis)
clim([physical_min_log, actual_max_log]); 
% ==============================================

axis equal;
xlabel('x'); ylabel('y');
shading interp; 
colorbar;