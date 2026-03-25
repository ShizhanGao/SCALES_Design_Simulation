% =========================================================================
% SCALES 阻力帆大修：动态仿真数据可视化 (Reviewer 2 - Major 4 & 7)
% 包含图A(Trade-off), 图B(微观时序), 图C(容错代偿)
% =========================================================================
clear; clc; close all;

% 全局字体与字号设置 (符合顶刊规范)
font_name = 'Times New Roman';
font_size = 12;
color_blue = '#0072BD'; % 高级蓝
color_red  = '#D95319'; % 砖红

%% ========================================================================
% 图 A (Fig. X1): 不同 n 情况下的 Trade-off 权衡折线图
% ========================================================================
figure('Name', 'Trade-off Analysis of Array Resolution', 'Position',[100, 100, 600, 400]);

% --- 【接口 1：你需要补充的真实数据】 ---
n_array =[100, 300, 400, 500];         % X轴：子帆数量
rmse_data =[0.085, 0.022, 0.015, 0.012]; % 左Y轴：稳态跟踪RMSE误差 (度)
freq_data =[0.2, 0.65, 1.1, 1.6];      % 右Y轴：单帆平均切换频率 (次/分钟)
% ----------------------------------------

yyaxis left
plot(n_array, rmse_data, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', color_blue);
ylabel('Steady-State Pointing RMSE (deg)', 'FontName', font_name, 'FontSize', font_size);
ylim([0, 0.1]); 
ax = gca; ax.YColor = color_blue; % 左轴颜色与数据对应

yyaxis right
plot(n_array, freq_data, '--s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', color_red);
ylabel('Average Switching Rate (flips/min/sail)', 'FontName', font_name, 'FontSize', font_size);
ylim([0, 2.5]);
ax = gca; ax.YColor = color_red; % 右轴颜色与数据对应

xlabel('Number of Discrete Elements (n)', 'FontName', font_name, 'FontSize', font_size);
title('Trade-off: Tracking Accuracy vs. Mechanical Wear', 'FontName', font_name, 'FontSize', font_size+2);
grid on; set(gca, 'FontName', font_name, 'FontSize', font_size);
legend({'Pointing RMSE', 'Switching Rate'}, 'Location', 'northeast', 'FontName', font_name);

%% ========================================================================
% 图 B (Fig. X2): n=500 时的典型子帆动作阶梯图 (Stairs Plot)
% ========================================================================
figure('Name', 'Microscopic Switching Behavior', 'Position',[150, 150, 600, 500]);

% --- 【接口 2：你需要补充的真实数据】 ---
time_vector = 0:1:140; % 假设仿真140秒，步长1s
% 假设从你的 Delta_hist 矩阵中挑出的 3 个代表性子帆状态向量 (只含 0 和 1)
% 伪造数据以作展示：
sail_heavy = randi([0,1], 1, length(time_vector)); % 频繁切换的
sail_mid   =[ones(1,40), zeros(1,20), ones(1,81)]; % 偶尔切换的
sail_lazy  =[zeros(1,100), ones(1,41)];            % 几乎不动的
% ----------------------------------------

% 子图 1：频繁微调的桨叶
subplot(3,1,1);
stairs(time_vector, sail_heavy, 'LineWidth', 1.5, 'Color', color_blue);
ylim([-0.2, 1.2]); yticks([0, 1]); xlim([0, max(time_vector)]);
ylabel('State', 'FontName', font_name);
title('Blade A: High-frequency Fine-tuning (Dynamic Load)', 'FontName', font_name);
set(gca, 'FontName', font_name, 'FontSize', 10); grid on;

% 子图 2：偶尔切换的桨叶
subplot(3,1,2);
stairs(time_vector, sail_mid, 'LineWidth', 1.5, 'Color', '#EDB120'); % 黄色
ylim([-0.2, 1.2]); yticks([0, 1]); xlim([0, max(time_vector)]);
ylabel('State', 'FontName', font_name);
title('Blade B: Occasional Transition (Posture Shift)', 'FontName', font_name);
set(gca, 'FontName', font_name, 'FontSize', 10); grid on;

% 子图 3：长期休眠的桨叶
subplot(3,1,3);
stairs(time_vector, sail_lazy, 'LineWidth', 1.5, 'Color', color_red);
ylim([-0.2, 1.2]); yticks([0, 1]); xlim([0, max(time_vector)]);
xlabel('Time (s)', 'FontName', font_name); ylabel('State', 'FontName', font_name);
title('Blade C: Prolonged Dormancy (Steady-state Anchoring)', 'FontName', font_name);
set(gca, 'FontName', font_name, 'FontSize', 10); grid on;

%% ========================================================================
% 图 C (Fig. X3): 故障率下的容错代偿机制 (双Y轴柱状图)
% ========================================================================
figure('Name', 'Fault-Tolerant Compensatory Mechanism', 'Position',[200, 200, 600, 400]);

% --- 【接口 3：你需要补充的真实数据】 ---
x_labels =[0, 10, 20, 30]; % 故障率 (%)
% 注意：这四个数据是你跑四次闭环容错仿真后提取的
rmse_fault  =[0.012, 0.013, 0.016, 0.045]; % 左Y轴：不同故障率下的稳态RMSE
freq_health = [1.6, 2.1, 3.8, 5.5];       % 右Y轴：*健康子帆*的平均切换频率
% ----------------------------------------

x_pos = 1:4;
bar_width = 0.35;

yyaxis left
bar(x_pos - bar_width/2, rmse_fault, bar_width, 'FaceColor', color_blue);
ylabel('Steady-State Pointing RMSE (deg)', 'FontName', font_name, 'FontSize', font_size);
ylim([0, 0.06]);
ax = gca; ax.YColor = color_blue;

yyaxis right
bar(x_pos + bar_width/2, freq_health, bar_width, 'FaceColor', color_red);
ylabel('Healthy Blade Switching Rate (flips/min/sail)', 'FontName', font_name, 'FontSize', font_size);
ylim([0, 7]);
ax = gca; ax.YColor = color_red;

xticks(x_pos);
xticklabels({'0% (Nominal)', '10% Fail', '20% Fail', '30% Fail'});
xlabel('Blade Failure Rate', 'FontName', font_name, 'FontSize', font_size);
title('Tolerance Evaluation: Accuracy vs. Compensatory Actuation', 'FontName', font_name, 'FontSize', font_size+2);
grid on; set(gca, 'FontName', font_name, 'FontSize', font_size);
legend({'Pointing RMSE', 'Healthy Blade Switching Rate'}, 'Location', 'northwest', 'FontName', font_name);