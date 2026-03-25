function Verify_Trajectory_ODE45(Z_ref, V_ref, I)


    %% 2. 参数设置
    dt = 0.4; % 步长
    N = size(V_ref, 2); % 控制步数
    
    % 如果 V_ref 比 Z_ref 少一列（常见情况），N 取 V_ref 的长度
    % 如果 Z_ref 比 V_ref 短，取 Z_ref 的长度
    if size(Z_ref, 2) < N + 1
        N = size(Z_ref, 2) - 1;
    end

    %% 3. ODE45 逐都积分 (Zero-Order Hold)
    % 我们不能一次性把 V_hist 给 ode45，必须一步一步来，
    % 因为 ode45 需要在 0.4s 内保持控制量恒定。
    
    X_ode = zeros(7, N+1);
    X_ode(:, 1) = Z_ref(:, 1); % 初始状态对齐
    
    x_curr = Z_ref(:, 1);
    
    fprintf('开始 ODE45 积分对比 (共 %d 步)...\n', N);
    
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9); % 高精度
    
    for k = 1:N
        % 当前时间段
        t_span = [(k-1)*dt, k*dt];
        
        % 当前控制量 (假设欠驱动，Z轴力矩为0)
        u_curr = V_ref(:, k); 
        u_full = [u_curr; 0]; 
        
        % 执行积分
        [~, x_temp] = ode45(@(t,y) dynamics(t, y, u_full, I), t_span, x_curr, options);
        
        % 取最后时刻状态作为下一步起点
        x_next = x_temp(end, :)';
        
        % 【关键】四元数归一化 (消除积分累积误差)
        x_next(1:4) = x_next(1:4) / norm(x_next(1:4));
        
        % 存储
        X_ode(:, k+1) = x_next;
        x_curr = x_next;
    end
    
    %% 4. 绘图对比
    t_axis = 0:dt:(N*dt);
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 700]);
    
    % --- 子图1: 四元数 ---
    subplot(2, 2, 1);
    hold on;
    plot(t_axis, Z_ref(1:4, 1:N+1), 'o', 'MarkerSize', 6, 'LineWidth', 1);
    set(gca, 'ColorOrderIndex', 1); % 重置颜色循环
    plot(t_axis, X_ode(1:4, :), '-', 'LineWidth', 1.5);
    title('四元数对比 (点: Z\_hist, 线: ODE45)');
    legend('q0_{ref}', 'q1_{ref}', 'q2_{ref}', 'q3_{ref}');
    grid on; xlabel('Time (s)'); ylabel('Quaternion');
    
    % --- 子图2: 角速度 ---
    subplot(2, 2, 2);
    hold on;
    plot(t_axis, Z_ref(5:7, 1:N+1), 'o', 'MarkerSize', 6, 'LineWidth', 1);
    set(gca, 'ColorOrderIndex', 1);
    plot(t_axis, X_ode(5:7, :), '-', 'LineWidth', 1.5);
    title('角速度对比 (点: Z\_hist, 线: ODE45)');
    legend('wx_{ref}', 'wy_{ref}', 'wz_{ref}');
    grid on; xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
    
    % --- 子图3: 误差分析 ---
    subplot(2, 2, [3, 4]);
    err_q = Z_ref(1:4, 1:N+1) - X_ode(1:4, :);
    err_w = Z_ref(5:7, 1:N+1) - X_ode(5:7, :);
    
    % 计算每一时刻的范数误差
    norm_err_q = sqrt(sum(err_q.^2, 1));
    norm_err_w = sqrt(sum(err_w.^2, 1));
    
    plot(t_axis, norm_err_q, 'b-o', 'LineWidth', 1.5); hold on;
    plot(t_axis, norm_err_w, 'r-^', 'LineWidth', 1.5);
    title('状态误差 (Z\_hist - ODE45)');
    legend('Quaternion Error Norm', 'Omega Error Norm');
    grid on; xlabel('Time (s)'); ylabel('Error Norm');
    
    % 打印最大误差
    fprintf('四元数最大误差: %.2e\n', max(norm_err_q));
    fprintf('角速度最大误差: %.2e\n', max(norm_err_w));
    
    if max(norm_err_q) < 1e-4
        fprintf('结论: 轨迹一致性良好。\n');
    else
        fprintf('结论: 存在较大偏差。可能原因：\n');
        fprintf('1. Z_hist 生成时使用的积分器精度较低 (如欧拉法或低阶RK)。\n');
        fprintf('2. 动力学模型参数 (I) 不匹配。\n');
    end
end

%% --- 动力学函数 ---
function dy = dynamics(~, y, u, I)
    q = y(1:4);
    w = y(5:7);
    
    % 运动学: dq = 0.5 * G(q) * w
    G = [-q(2) -q(3) -q(4);
          q(1) -q(4)  q(3);
          q(4)  q(1) -q(2);
         -q(3)  q(2)  q(1)];
    dq = 0.5 * G * w;
    
    % 动力学: I * dw = u - w x (Iw)
    % dw = I \ (u - cross(w, I*w));
    Iw = I * w;
    dw = I \ (u - cross(w, Iw));
    
    dy = [dq; dw];
end

%% --- 生成测试数据 (仅当缺数据时使用) ---
function [Z, V, I] = generate_dummy_data()
    I = diag([10; 20; 30]);
    dt = 0.4;
    N = 50;
    V = 0.1 * sin(0.1 * (1:N)); % 简单的正弦控制
    V = [V; V]; % 2维控制
    
    Z = zeros(7, N+1);
    Z(:,1) = [1; 0; 0; 0; 0.1; 0.1; 0.1];
    
    % 用简单的 RK4 生成 Z_hist
    x = Z(:,1);
    for k = 1:N
        u = [V(:,k); 0];
        % RK4 Step
        f = @(xx, uu) dynamics(0, xx, uu, I);
        k1 = f(x, u);
        k2 = f(x + 0.5*dt*k1, u);
        k3 = f(x + 0.5*dt*k2, u);
        k4 = f(x + dt*k3, u);
        x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        x(1:4) = x(1:4)/norm(x(1:4));
        Z(:, k+1) = x;
    end
end