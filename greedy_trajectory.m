% att_mpc_casadi_debug_fixed.m
% CasADi NMPC for underactuated attitude (quaternion)

clear; clc; close all;
import casadi.*
load('para.mat')

sc.n  = 200;        % n个子帆 (较大数量更能看出包络形状)
params_Code;

%% ---------- ref. traj. parameters ----------
nu = 2; nx = 7; % 控制量与状态量个数

Ir = sc.J; % 转动惯量矩阵
% Ir = diag(Ir([2, 3, 1]));
% I = diag([0.13; 0.1198; 0.0144]); 
Iinv = inv(Ir);

% 
% u_max = 0.06; u_min = -0.06;
u_set = 0.06;
u_max = 0.8*u_set; u_min = -0.8*u_set;

% 采样间隔
dtr = 3;
% 预测步长
Nr = 40;

% 初始状态参数
q0 = eul2quat(deg2rad([40 25 10]), 'ZYX')';
% q0 = quatnormalize([1 0.17 0 -0.02])';
w0 = deg2rad([0.01; -0.05; 0.05]);


x0 = [q0; w0];
x_des = [1;0;0;0; 0; 1.1455e-3; 0];

% 二次型参数 
Wq = 1; Ww = 10; Wu = 1;

% 终端状态
Zf = [0.05,0.005,0.05,0.005,0.0005,0.005]' * 1e-3;
lbg = [0.995; -Zf];
ubg = [1; Zf];
I3 = SX.eye(3);
%% 工况设置
It = diag([1.1 0.9 1.05])*Ir;
% It = Ir;



%% ---------- CasADi 符号计算设置 ----------
u = SX.sym('u', nu*Nr);    % decision vector: [u0x u0y u1x u1y ...]'

% 转换转动惯量矩阵至符号计算
I_c = SX(Ir); Iinv_c = SX(Iinv);

% 初始状态的符号表示，向量(7维)
x_sym = SX.sym('x0', nx);

cost = 0;
xk = x_sym;

for i = 1:Nr % 预测步长
    ui = u((i-1)*nu+1:(i-1)*nu+nu);
    ufull = [SX(0); ui];  % underactuated: z torque zero
    
    % RK4 single step using CasADi operations
    xk = casadi_rk4_sym(xk, ufull, dtr, I_c, Iinv_c);
    % normalize quaternion symbolically
    qk = xk(1:4);
    wk = xk(5:7);
    qk = qk / sqrt(sumsqr(qk));
    xk = [qk; wk];
    
    % quaternion error qe = q_des^* * qk
    qe = quat_mul_cas(quat_conj_cas(x_des(1:4)), qk);

    e0 = qe(1); 
    e_vec = qe(2:4);
    
    % (B) 角速度误差 (带旋转矩阵修正)
    skew_e = [0, -e_vec(3), e_vec(2); e_vec(3), 0, -e_vec(1); -e_vec(2), e_vec(1), 0];
    R_matrix = (e0^2 - dot(e_vec, e_vec)) * I3 + 2 * (e_vec * e_vec') - 2 * e0 * skew_e;
    err_w = wk - mtimes(R_matrix, x_des(5:7));

    % cost terms
    err_q = (1 - e0)^2;
    cost = cost + Wq * err_q + Ww * (err_w' * err_w) + Wu * (ui' * ui);

end
constraint_term = [qe; err_w];
nlp = struct('x', u, 'f', cost, 'g', constraint_term, 'p', x_sym);
% constraint_term = Wq * err_q + Ww * (err_w' * err_w);
% Build NLP (u is decision, x_sym is parameter)
% nlp = struct('x', u, 'f', cost, 'g', constraint_term, 'p', x_sym);
% nlp = struct('x', u, 'f', cost, 'p', x_sym);

% IPOPT options (moderate verbosity)
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = false;
opts.ipopt.max_iter = 500;
opts.ipopt.tol = 1e-9;
solver = nlpsol('solver','ipopt', nlp, opts);

%% ---------- setup for closed-loop ----------
lb = u_min * ones(nu*Nr,1);
ub = u_max * ones(nu*Nr,1);
u0 = zeros(nu*Nr,1);

N_add = 5;
Zhist = zeros(nx, N_add+Nr+1); Zhist(:,1) = x0;
Vhist = zeros(nu, N_add+Nr);
x = x0;

fprintf("\n===== CASADI NMPC (fixed quaternion ops) START =====\n\n");

tic
sol = solver('x0', u0, 'lbx', lb, 'ubx', ub, 'lbg', lbg, 'ubg', ubg, 'p', x);
fprintf(" iteration: %d\n", i);
toc
% --- 新增检查代码 ---
stats = solver.stats();
if ~strcmp(stats.return_status, 'Solve_Succeeded')
    warning('MPC求解失败! 状态: %s', stats.return_status);
    % 此时 sol.x 里的结果是不可信的，通常是不满足约束的
end

u_opt = full(sol.x);
fval = full(sol.f);
stats = solver.stats();

% solver.stats may be nested struct or struct array, handle common fields
if isfield(stats,'return_status'); status = stats.return_status; else status = 'unknown'; end
if isfield(stats,'iter_count'); iters = stats.iter_count; else iters = NaN; end

fprintf(" solver status: %s\n", string(status));
fprintf(" iterations: %d\n", iters);
fprintf(" final cost: %.6e\n", fval);
fprintf(" optimal sequence norm = %.3e\n", norm(u_opt));

% check bounds touch
if any(abs(u_opt) >= 0.999*u_max)
    fprintf(" NOTE: some controls in optimal sequence touch bounds.\n");
end

% end
fprintf("\n");

fprintf("===== MPC finished =====\n");

for k = 1:(Nr+N_add)
        if k > Nr
            u_apply = [0 0]';
        else
            u_apply = u_opt((k-1)*nu+1: k*nu);
        end
        
        % simulate one step on true plant (RK4 numeric)
        x = rk4_step_numeric(x, [0; u_apply], dtr, Ir, Iinv);
        % normalize numerically quaternion
        x(1:4) = x(1:4) / norm(x(1:4));
        
        % store
        Zhist(:,k+1) = x;
        Vhist(:,k) = u_apply;
        % Jhist(:,k) = fval;
end
%% ---------- plots ----------
time = dtr*(0:Nr+N_add);

quat_to_euler = @(q) quat2eul(q'); % needs q as row [q0 q1 q2 q3], returns [yaw pitch roll]
euler_hist = zeros(3, Nr+1);
for i=1:Nr+N_add+1
    eul = quat_to_euler(Zhist(1:4,i)); % returns ZYX: [yaw pitch roll]
    euler_hist(:,i) = eul';
end
% plot(time, rad2deg(euler_hist(1,:)),'-','LineWidth',1.2); hold on;
% plot(time, rad2deg(euler_hist(2,:)),'-','LineWidth',1.2);
% plot(time, rad2deg(euler_hist(3,:)),'-','LineWidth',1.2);
% legend('yaw (z)','pitch (y)','roll (x)'); xlabel('t (s)'); ylabel('deg'); grid on;

%% ---------- Ancillary Controller ----------
% It = Ir; % Ir*diag([0.95 1.05 1.01]); % 真实转动惯量矩阵
Iinv = inv(It);

% u_max = 0.06; u_min = -0.06;
u_max = u_set; u_min = -u_set;

% 采样间隔
dt = 1;
% 预测步长
N = 20;

% 终端状态
Eps_f = [2 2 5 0.1 0.1 0.05]'*1e-3;
lbg = [0.995; -Eps_f];
ubg = [1; Eps_f];

% ---------- CasADi 符号计算设置 (辅助控制器) ----------
solver_anc = build_ancillary_mpc(Ir, dt, N);

%% --- 构建 NLP 仿真参数设置 ---
% 注意 p 现在包含了初始状态 + 参考轨迹序列

TotalTime = (Nr + N_add) * dtr;   % 总时间
t_ref_arr = dtr * (0:Nr + N_add); % 参考轨迹时间轴

% 仿真总步数
SimSteps = floor(TotalTime / dt);

% 初始化记录数组
log_X = zeros(nx, SimSteps+1);
log_Ua = zeros(nu, SimSteps);
log_U  = zeros(nu+1, SimSteps);
log_Cost = zeros(1, SimSteps);
log_Time = 0:dt:TotalTime;

% 初始状态 (加入一点点初始偏差以测试稳定性)
x_real = Zhist(:,1); 
log_X(:, 1) = x_real;

% Warm Start 初值
u_last = [];

fprintf('==============================================\n');
fprintf('开始 Tube-based MPC 跟踪仿真\n');
fprintf('总时长: %.2fs, 仿真步长: %.2fs, 总步数: %d\n', TotalTime, dt, SimSteps);
fprintf('==============================================\n');

%% --- 循环开始 ---

    % ----------------- 添加了操纵律 -----------------
    Delta_init =[zeros(sc.n/2, 1); ones(sc.n/2, 1)];
    q_err = [1 0 0 0]';
    M_array = Calculate_SCALES_plot2(sc, env, q_err, Delta_init);
    
    % 提取 My 和 Mz 并计算绝对值参数
    My = M_array(:, 2);
    Mz = M_array(:, 3);
    
    c_val = abs(My(1));
    a_val = abs(Mz(2) - Mz(1));
    b_val = abs(Mz(1)) - a_val;
    n = sc.n;
    m = n/4;
    % ------------------------------------------------
    p_fail = 0;
    sail_domains = cell(4*m, 1);
    for idx = 1:4*m
        sail_domains{idx} = get_sub_sail_domain(p_fail);
    end
    
    % ========================================================
    % 步骤 B: 提取宏观系流域 (根据最新映射规则推导)
    % ========================================================
    domain_u = cell(m, 1);
    domain_v = cell(m, 1);
    for i = 1:m
        % u 的映射: i 和 375+i (映射逻辑: u = s1 - s2)
        dom_u1 = sail_domains{i};
        dom_u2 = sail_domains{3*m + i};
        u_opts =[];
        for s1 = dom_u1, for s2 = dom_u2, u_opts(end+1) = s1 - s2; end, end %#ok<*AGROW>
        domain_u{i} = unique(u_opts);
        
        % v 的映射: 125+i 和 250+i (映射逻辑: v=1对应0 1, 即 v = s2 - s1)
        dom_v1 = sail_domains{m + i};
        dom_v2 = sail_domains{2*m + i};
        v_opts =[];
        for s1 = dom_v1, for s2 = dom_v2, v_opts(end+1) = s2 - s1; end, end % 注意这里的 s2 - s1
        domain_v{i} = unique(v_opts);
    end

Delta_hist = [];


for k = 1:SimSteps
    t_curr = (k-1) * dt;

    % 1. 构建参数 p (包含 x0 和未来 N+1 个参考点)
    % ref_seq 结构: [z_ref_0; u_ref_0; z_ref_1; u_ref_1; ... ; z_ref_N; u_ref_N]
    ref_seq_val = [];

    for j = 0:N
        t_pred = t_curr + j * dt;
        [z_ref_j, u_ref_j] = interpolate_reference(t_pred, t_ref_arr, Zhist, Vhist);
        ref_seq_val = [ref_seq_val; z_ref_j; u_ref_j];
        if k == 1
            u_last = [u_last; u_ref_j];
        end
    end
    u_last = u_last(1:nu*N);
    p_val = [x_real; ref_seq_val];

    % 2. 设置 Bounds (推力约束)
    % u_max 是标量或向量，假设为 0.06
    args.lbx = repmat([-u_max; -u_max], N, 1);
    args.ubx = repmat([ u_max;  u_max], N, 1);

    % 末端约束 lbg/ubg (必须与你构建solver时的维度一致)
    % 假设你之前的 lbg 是 7维 [0.995; -Eps_f]
    args.lbg = lbg; 
    args.ubg = ubg;

    % 3. 初始猜测 (Warm Start)
    args.x0 = u_last; 
    args.p  = p_val;
    
    % 4. 调用求解器
    sol = solver_anc('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, ...
                     'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
    fprintf(" iteration: %d\n", solver_anc.stats().iter_count);
    % 5. 检查求解状态
    stats = solver_anc.stats();
    if ~strcmp(stats.return_status, 'Solve_Succeeded')
        % 简单处理：如果求解失败，尝试使用参考控制量作为保底
        warning('Step %d (%.2fs): Infeasible! Applying Feedforward U.', k, t_curr);
        % u_opt_seq = ref_seq_val(nx+1:nx+nu); % 取参考控制
        u_opt_seq = full(sol.x);
        u_last = u_opt_seq;
        log_Cost(k) = full(sol.f);
    else
        u_opt_seq = full(sol.x);
        u_last = u_opt_seq; % 更新 warm start
        log_Cost(k) = full(sol.f);
    end
    % u_opt_seq = ref_seq_val(nx+1:nx+nu); % 取参考控制
    % 提取第一步控制量
    u_apply = u_opt_seq(1:nu);
    log_Ua(:, k) = u_apply;

    % 6. 真实环境演化 (Real Plant Simulation)
    % 使用 ode45 或 rk4_numeric，注意这里要用真实惯量 It 和干扰
    % 模拟扰动 (例如 1e-4 Nm 的常值+随机扰动)
    dist_torque = 0*randn(3,1); 
    
    % 积分一步
    % u_act = add_random_deflection([0; u_apply], 3);

    
    
    [u_calc, v_calc, ~, ~] = SCALES_FaultTolerant_Solver(a_val, b_val, c_val, u_apply(1), u_apply(2), domain_u, domain_v);
    
    % (4) 生成 500 位严格受限的物理控制指令
    Delta = Generate_FaultAware_Delta(u_calc, v_calc, sail_domains);
    Delta_hist = [Delta_hist; Delta'];
    % (5) 真实物理引擎解算
    M_out = Calculate_SCALES_Torque(sc, env, q_err, Delta);
    T_My(k) = M_out(2);
    T_Mz(k) = M_out(3);
    disp(norm(M_out(2:3) - u_apply))

    u_act = M_out;

    % u_act = add_random_deflection(M_out, 3);
    u_act = (0.8 + 0.4*rand(1))*u_act;
    [~, x_next_ode] = ode45(@(t,y) dynamics_numeric_with_dist(t, y, u_act, It, dist_torque), ...
                            [t_curr, t_curr + dt], x_real);
    log_U(:, k) = u_act;
    x_real = x_next_ode(end, :)';
    x_real(1:4) = x_real(1:4) / norm(x_real(1:4)); % 物理归一化
    log_X(:, k+1) = x_real;

    % 7. 打印进度
    if mod(k, 10) == 0 || k == 1
        % 计算当前跟踪误差范数
        [z_ref_curr, ~] = interpolate_reference(t_curr, t_ref_arr, Zhist, Vhist);
        err_norm = norm(x_real - z_ref_curr);

        fprintf('Step %3d/%d | T: %.2fs | Cost: %.2e | TrackErr: %.2e | Status: %s\n', ...
            k, SimSteps, t_curr, log_Cost(k), err_norm, stats.return_status);
    end
end
fprintf('仿真结束。\n');

%% --- 绘图前的数据预处理 ---
% 1. 将四元数转为欧拉角 (ZYX顺序: Yaw->Pitch->Roll)
% 注意: quat2eul 返回的是 [Yaw, Pitch, Roll]
eul_real = quat2eul(log_X(1:4,:)', 'ZYX') * 180/pi; % [deg]
eul_nom  = quat2eul(Zhist(1:4,:)', 'ZYX') * 180/pi; % [deg]

% 2. 对齐时间轴
% 假设 Zhist 是按 dtr 采样的离散点，log_Time 是仿真细密点
t_nom = 0:dtr:(size(Zhist,2)-1)*dtr; 

%% --- 仿真结果图 ---

fig3 = figure('Units', 'centimeters', 'Position', [5, 5, 22, 16], 'Color', 'w');
t1 = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 配色方案
c_real = [0.1, 0.3, 0.7];       % 实线颜色 (深蓝)
c_nom  = [0.85, 0.33, 0.1];     % 标称轨迹颜色 (深橘红)
c_unact = [0.5, 0.2, 0.6];      % 欠驱动轴特殊颜色 (紫色)

% (a) Euler Angles (三轴姿态)
ax1 = nexttile([1, 2]); % 跨两列
hold on;
% 绘制实际轨迹
plot(time, Zhist(1,:)', '--', 'Color', 'k', 'LineWidth', 1.5, 'DisplayName', '$q_{e0}$');
plot(time, Zhist(2,:)', 'Color', c_unact, 'LineWidth', 1.5, 'DisplayName', '$q_{e1}$');
plot(time, Zhist(3,:)', '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5, 'DisplayName', '$q_{e2}$');
plot(time, Zhist(4,:)', '--', 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', '$q_{e3}$');
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('quaternion', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Box', 'off', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');

ax1 = nexttile([1, 2]); % 跨两列
hold on;
% 绘制实际轨迹
plot(time, rad2deg(euler_hist(3,:)), 'Color', c_unact, 'LineWidth', 1.5, 'DisplayName', 'Roll (\phi) [Underactuated]');
plot(time, rad2deg(euler_hist(2,:)), '--', 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5, 'DisplayName', 'Pitch (\theta)');
plot(time, rad2deg(euler_hist(1,:)), '--', 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', 'Yaw (\psi)');

ylabel('Angle (deg)', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Box', 'off', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
% title('(a) Attitude Stabilization (Solid: Real, Dashed: Nominal)', 'FontName', 'Times New Roman');

% (b) Actuated Angular Velocities (Y & Z axes)
ax2 = nexttile([1, 2]);
hold on;
% Y轴 (Pitch)
plot(time, Zhist(5,:)', 'Color', c_unact, 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_x$');
plot(time, Zhist(6,:)', 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_y$');
% Z轴 (Yaw)
plot(time, Zhist(7,:)', 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_z$');

ylabel('$\omega$ (rad/s)', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
% title('(b) Angular Velocity (Actuated Axes)', 'FontName', 'Times New Roman');

% (c) Underactuated Angular Velocity (X axis) - 重点关注对象
ax3 = nexttile([1, 2]);
hold on;
stairs(time, [Vhist(1,:),Vhist(1,end)]', 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.5, 'DisplayName', '$u_y$'); 
grid on
stairs(time, [Vhist(2,:),Vhist(2,end)]', 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', '$u_z$');
% Z轴 (Yaw)

ylabel('Torque $u$ (Nm)', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
% title('(c) Angular Velocity (Underactuated Axis)', 'FontName', 'Times New Roman');

% --- 图 1: 姿态与角速度 (Attitude & Angular Velocity) ---
fig1 = figure('Units', 'centimeters', 'Position', [5, 5, 22, 16], 'Color', 'w');
t1 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 配色方案
c_real = [0.1, 0.3, 0.7];       % 实线颜色 (深蓝)
c_nom  = [0.85, 0.33, 0.1];     % 标称轨迹颜色 (深橘红)
c_unact = [0.5, 0.2, 0.6];      % 欠驱动轴特殊颜色 (紫色)

% (a) Euler Angles (三轴姿态)
ax1 = nexttile([1, 2]); % 跨两列
hold on;
% 绘制实际轨迹
plot(log_Time, eul_real(:,3), 'Color', c_unact, 'LineWidth', 1.5, 'DisplayName', 'Roll (\phi) [Underactuated]');
plot(log_Time, eul_real(:,2), 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.5, 'DisplayName', 'Pitch (\theta)');
plot(log_Time, eul_real(:,1), 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', 'Yaw (\psi)');
% 绘制标称轨迹 (点划线)
plot(t_nom, eul_nom(:,3), '--', 'Color', c_unact, 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(t_nom, eul_nom(:,2), '--', 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot(t_nom, eul_nom(:,1), '--', 'Color', c_real, 'LineWidth', 1.2, 'HandleVisibility', 'off');

ylabel('Angle (deg)', 'Interpreter', 'latex');
legend('Location', 'northeast', 'Box', 'off');
grid on; set(gca, 'FontName', 'Times New Roman');
title('(a) Attitude Stabilization (Solid: Real, Dashed: Nominal)', 'FontName', 'Times New Roman');

% (b) Actuated Angular Velocities (Y & Z axes)
ax2 = nexttile;
hold on;
% Y轴 (Pitch)
plot(log_Time, log_X(6,:), 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.5, 'DisplayName', 'Real $\omega_y$');
plot(t_nom, Zhist(6,:), '--', 'Color', c_nom, 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_y$');
% Z轴 (Yaw)
plot(log_Time, log_X(7,:), 'Color', c_real, 'LineWidth', 1.5, 'DisplayName', 'Real $\omega_z$');
plot(t_nom, Zhist(7,:), '--', 'Color', c_nom, 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_z$');

ylabel('$\omega_{y,z}$ (rad/s)', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
title('(b) Angular Velocity (Actuated Axes)', 'FontName', 'Times New Roman');

% (c) Underactuated Angular Velocity (X axis) - 重点关注对象
ax3 = nexttile;
hold on;
plot(log_Time, log_X(5,:), 'Color', c_unact, 'LineWidth', 1.5, 'DisplayName', 'Real $\omega_x$');
plot(t_nom, Zhist(5,:), '--', 'Color', c_nom, 'LineWidth', 1.5, 'DisplayName', 'Nominal $\omega_x$');

ylabel('$\omega_x$ (rad/s)', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
title('(c) Angular Velocity (Underactuated Axis)', 'FontName', 'Times New Roman');


% --- 图 2: 控制力矩与误差 (Control Inputs & Tube Error) ---
fig2 = figure('Units', 'centimeters', 'Position', [10, 5, 22, 14], 'Color', 'w');
t2 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Control Torque u_y (Pitch) - 对应 log_U 的第1行
ax4 = nexttile([1, 2]);
hold on;
yline(u_max, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
yline(u_min, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
% 实际控制
stairs(log_Time(1:end-1), log_U(1,:), 'Color', c_unact, 'LineWidth', 1.2, 'DisplayName', '$u_x$');
stairs(log_Time(1:end-1), log_U(2,:), 'Color', [0.2, 0.6, 0.3], 'LineWidth', 1.2, 'DisplayName', '$u_y$');
% 标称控制
stairs(log_Time(1:end-1), log_Ua(1,:), '--', 'Color', c_nom, 'LineWidth', 1.2, 'DisplayName', 'Nominal $\bar{u}_y$');

ylim([-1.2*u_max, 1.2*u_max])
ylabel('Torque $u_y$ (Nm)', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
title('(a) Control Effort: Pitch Axis', 'FontName', 'Times New Roman');
hold on;
stairs(log_Time(1:end-1), log_U(3,:), 'Color', c_real, 'LineWidth', 1.2, 'DisplayName', '$u_z$');
stairs(log_Time(1:end-1), log_Ua(2,:), '--', 'Color', c_nom, 'LineWidth', 1.2, 'DisplayName', 'Nominal $\bar{u}_z$');

ylabel('Torque $u$ (Nm)', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
legend('Location', 'best', 'Box', 'off', 'Interpreter', 'latex');
title('(a) ', 'FontName', 'Times New Roman');

% (c) Tube Error Norm (误差范数)
ax6 = nexttile([1, 2]);
hold on;
err_norm = zeros(1, length(log_Time));
for i = 1:length(log_Time)
    % 简单的最近邻插值对齐时间
    [~, idx] = min(abs(log_Time(i) - t_nom));
    % 计算误差 (排除插值边缘效应)
    if idx <= size(Zhist, 2)
        z_curr = Zhist(:, idx);
        err_norm(i) = norm(log_X(2:4, i) - z_curr(2:4));
    end
end
area(log_Time, err_norm, 'FaceColor', [0.85, 0.9, 0.95], 'EdgeColor', [0.2, 0.5, 0.7], 'LineWidth', 1);

ylabel('Error $\|\mathbf{x} - \bar{\mathbf{x}}\|$', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
grid on; set(gca, 'FontName', 'Times New Roman');
title('(b)', 'FontName', 'Times New Roman');

% 标注说明
% Add annotation to explain Underactuation
% annotation(fig2, 'textbox', [0.75, 0.85, 0.2, 0.05], 'String', 'Control u_x \equiv 0', ...
%     'EdgeColor', 'none', 'FontName', 'Times New Roman', 'FontSize', 10, 'Interpreter', 'tex');