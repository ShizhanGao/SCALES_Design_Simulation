% 测试脚本
clc; clear; close all;

v0 = [-1; 0; 0]; % 假设这是你的气流方向
angle_err = 3;   % 3度偏差

% 1. 单次测试
v_new = add_random_deflection(v0, angle_err);
actual_angle = rad2deg(acos(dot(v0, v_new) / (norm(v0)*norm(v_new))));
fprintf('设定偏差: %.2f°, 实际计算偏差: %.4f°\n', angle_err, actual_angle);

% 2. 蒙特卡洛可视化 (画个圈看看是不是圆锥)
N = 500;
V_cloud = zeros(3, N);
for i = 1:N
    V_cloud(:, i) = add_random_deflection(v0, angle_err);
end

figure;
quiver3(0,0,0, v0(1), v0(2), v0(3), 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5); hold on;
scatter3(V_cloud(1,:), V_cloud(2,:), V_cloud(3,:), 10, 'b', 'filled');
axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Random 3^\circ Deflection Cone');
legend('Original Vector', 'Deflected Samples');
view(3);