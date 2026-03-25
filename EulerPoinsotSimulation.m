function EulerPoinsotSimulation()
    % 欧拉-潘索运动可视化 (Euler-Poinsot Motion Visualization)
    % 展示刚体坐标系下运动的闭合性 vs 惯性坐标系下运动的非闭合性
    
    clc; clear; close all;

    %% 1. 系统参数设置
    % 选取不对称刚体 (I1 < I2 < I3)
    I = [10; 20; 40]; 
    I_mat = diag(I);
    I_inv = inv(I_mat);

    % 初始条件
    % 初始角速度 (选取一个一般位置，避开主轴，以产生明显的进动/章动)
    w0 = [0.8; 0.8; 8]; 
    
    % 初始姿态 (四元数 [qw, qx, qy, qz])
    q0 = [1; 0; 0; 0];
    
    % 仿真时间 (设置较长时间以观察非闭合性)
    T_span = [0, 1500]; 
    
    %% 2. 动力学积分 (ODE45)
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    init_state = [q0; w0];
    [t, state] = ode45(@(t,y) rigid_body_dynamics(t, y, I, I_inv), T_span, init_state, options);

    % 提取结果
    Q = state(:, 1:4);   % 四元数
    W = state(:, 5:7);   % 刚体系角速度
    
    %% 3. 数据处理
    
    % --- 计算能量和角动量 (用于绘制参考曲面) ---
    % 动能 T = 0.5 * w' * I * w
    E_k = 0.5 * (I(1)*w0(1)^2 + I(2)*w0(2)^2 + I(3)*w0(3)^2);
    % 角动量大小 L^2 = |I*w|^2
    L_vec_0 = I_mat * w0;
    L_mag = norm(L_vec_0);
    
    % --- 计算惯性系下的轨迹 ---
    % 我们追踪刚体的 Z 轴在惯性系单位球上的投影
    BodyZ_axis = [0; 0; 1];
    Inertial_Trace = zeros(length(t), 3);
    
    for k = 1:length(t)
        % 将四元数转换为旋转矩阵 R_IB
        q_curr = Q(k, :);
        % 这里使用简化的四元数转旋转矩阵逻辑 (或者使用 quat2rotm 如果有工具箱)
        R = quat_to_rot(q_curr);
        
        % 将刚体 Z 轴转到惯性系: v_inertial = R * v_body
        Inertial_Trace(k, :) = (R * BodyZ_axis)';
    end

    %% 4. 绘图可视化
    
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    
    % === 子图 1: 刚体坐标系 (Body Frame) ===
    subplot(1, 2, 1);
    hold on; grid on; axis equal;
    xlabel('\omega_x'); ylabel('\omega_y'); zlabel('\omega_z');
    title({'【刚体坐标系】角速度轨迹', '结论: 闭合曲线 (Polhode)'}, 'FontSize', 12, 'FontWeight', 'bold');
    
    % 绘制角动量球 (L = |Iw|) => (Ix*wx)^2 + ... = L^2
    % 注意：这实际上是一个关于 w 的椭球，但在角动量空间是球。
    % 让我们直接画 w 所在的两个椭球面的交线原理
    
    % 1. 动能椭球 (Energy Ellipsoid): Ix*wx^2 + Iy*wy^2 + Iz*wz^2 = 2E
    [sx, sy, sz] = sphere(50);
    ex = sx * sqrt(2*E_k/I(1));
    ey = sy * sqrt(2*E_k/I(2));
    ez = sz * sqrt(2*E_k/I(3));
    surf(ex, ey, ez, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'b');
    
    % 2. 角动量守恒面 (Intersection surface for L conservation)
    % (Ix*wx)^2 + (Iy*wy)^2 + (Iz*wz)^2 = L^2
    lx = sx * (L_mag/I(1)); % 这不是正球体，这又是一个椭球
    ly = sy * (L_mag/I(2));
    lz = sz * (L_mag/I(3));
    % 为了视觉清晰，我们只画轨迹，用文字说明原理
    
    % 绘制 w 的轨迹
    plot3(W(:,1), W(:,2), W(:,3), 'r', 'LineWidth', 2);
    plot3(W(1,1), W(1,2), W(1,3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k'); % 起点
    
    view(135, 30);
    legend('动能椭球', '角速度轨迹', '起点');
    
    % === 子图 2: 惯性坐标系 (Inertial Frame) ===
    subplot(1, 2, 2);
    hold on; grid on; axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title({'【惯性坐标系】刚体Z轴指向', '结论: 准周期/不闭合 (Herpolhode-like)'}, 'FontSize', 12, 'FontWeight', 'bold');
    
    % 绘制单位球参考面
    surf(sx, sy, sz, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);
    
    % 绘制轨迹
    % 使用渐变色表示时间流逝，更能看出"不闭合"的缠绕感
    patch([Inertial_Trace(:,1); nan], [Inertial_Trace(:,2); nan], [Inertial_Trace(:,3); nan], ...
          [t; nan], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 1.5);
    colormap(jet);
    cb = colorbar; cb.Label.String = 'Time (s)';
    
    plot3(Inertial_Trace(1,1), Inertial_Trace(1,2), Inertial_Trace(1,3), 'ko', 'MarkerFaceColor', 'k');
    
    view(135, 30);
end

%% --- 辅助函数 ---

function dY = rigid_body_dynamics(~, Y, I, I_inv)
    % Y = [q0; q1; q2; q3; wx; wy; wz] (7x1)
    q = Y(1:4);
    w = Y(5:7);
    
    % 1. 归一化四元数 (防止数值漂移)
    q = q / norm(q);
    
    % 2. 欧拉动力学方程 (Body Frame)
    % I * dw/dt + w x (I*w) = 0  =>  dw/dt = -I_inv * (w x (I*w))
    Iw = I .* w; % 因为 I 是对角阵，可以直接点乘
    dw = -I_inv * cross(w, Iw);
    
    % 3. 运动学方程 (Quaternions)
    % dq/dt = 0.5 * q * [0; w]
    % 构造四元数乘法矩阵
    w_vec = [0; w];
    dq = 0.5 * quat_mul(q, w_vec);
    
    dY = [dq; dw];
end

function r = quat_mul(p, q)
    % 四元数乘法 r = p * q
    p0 = p(1); pv = p(2:4);
    q0 = q(1); qv = q(2:4);
    
    r0 = p0*q0 - dot(pv, qv);
    rv = p0*qv + q0*pv + cross(pv, qv);
    r = [r0; rv];
end

function R = quat_to_rot(q)
    % 四元数转旋转矩阵 R_IB (Body to Inertial)
    % q = [w, x, y, z]
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    R = [1 - 2*y^2 - 2*z^2,   2*x*y - 2*z*w,     2*x*z + 2*y*w;
         2*x*y + 2*z*w,       1 - 2*x^2 - 2*z^2, 2*y*z - 2*x*w;
         2*x*z - 2*y*w,       2*y*z + 2*x*w,     1 - 2*x^2 - 2*y^2];
end