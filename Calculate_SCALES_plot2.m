function M_array = Calculate_SCALES_plot2(para, env, q_err, Delta)
% CALCULATE_SCALES_PLOT2 计算分布式大气帆的气动合力矩与合力
%
% 输入:
%   para         : 结构体, 包含几何参数 (n, s, lb, wb, m1...m4)
%                  必须与 Calculate_SCALES_Inertia 的输入一致
%   env          : 结构体, 包含环境参数
%                  .h (轨道高度 m), .rho (大气密度), .Tw, .Tinf, .M (摩尔质量), .sigma_n, .sigma_t
%   q_err        : [4x1] 误差四元数 [q0; q1; q2; q3] (标量在前), 描述 S_o 到 S_b 的旋转
%   delta_angles : [nx1] 各子帆的转角 (弧度), 0表示法线沿x轴
%
% 输出:
%   M_b          : [3x1] 本体系下的总气动力矩 [N*m]
%   F_b          : [3x1] 本体系下的总气动力 [N]

    %% 1. 基础参数提取与几何生成
    n = para.n;
    s = para.s;
    wb = para.wb;
    lb = para.lb;
    
    % 计算单个子帆面积
    Area_i = lb * wb;
    
    % --- 生成所有子帆的位置矢量 r_i ---
    % 假设子帆分布在四个象限 (与惯量计算逻辑一致)
    % 布局：上下两行 (z = +/- lb/2)，左右对称分布
    % y 坐标: s/2 + (k-0.5)*wb
    
    r_pos = zeros(3, n); % 存储位置 [x; y; z]
    cols_per_quad = n / 4;
    
    idx = 1;
    % 遍历四个象限生成坐标
    for row_sign = [1, -1] % 下(-1), 上(+1)
        z_c = row_sign * (lb / 2); 
        for side_sign = [-1, 1] % 左(-1), 右(+1)
            for k = 1:cols_per_quad
                y_c = side_sign * (s/2 + (k - 0.5) * wb);
                x_c = 0; % 假设帆面平整，位于 x=0 平面
                
                r_pos(:, idx) = [x_c; y_c; z_c];
                idx = idx + 1;
            end
        end
    end
    
    %% 2. 速度矢量转换 (Orbital -> Body)
    % 计算圆轨道速度
    Re = 6371e3; % 地球半径
    mu = 3.986004418e14;
    r_orbit = Re + env.h;
    V_mag = sqrt(mu / r_orbit);
    
    % 轨道系下的来流速度 (假设相对速度即为负的轨道速度)
    % 轨道系定义: x_o 沿速度方向。故来流为 [-V, 0, 0]
    v_o = [-V_mag; 0; 0];
    
    % 利用误差四元数构建旋转矩阵 C_bo (Orbital to Body)
    % q_err = [q0, qv]'
    q0 = q_err(1);
    qv = q_err(2:4);
    
    % Rodrigues 公式 / 四元数转DCM
    I3 = eye(3);
    qv_cross = [0, -qv(3), qv(2); qv(3), 0, -qv(1); -qv(2), qv(1), 0];
    C_bo = (q0^2 - qv'*qv)*I3 + 2*(qv*qv') - 2*q0*qv_cross;
    
    % 本体系下的相对速度
    v_b = C_bo * v_o;
    
    % 速度单位矢量与大小
    v_rel_mag = norm(v_b);
    u_flow = v_b / v_rel_mag; % 本体系下来流方向单位矢量
    
    %% 3. 循环计算每个子帆的受力
    F_hist = zeros(n,3);
    M_hist = zeros(n,3);

    for i = 1:n
        if Delta(i) == 0
            dni = [1 0 0]';
        else
            dni = [0 1 0]';
        end
        
        ni = sgn(dot(u_flow, dni))*dni;

        v_tan = u_flow - dot(u_flow, ni)*ni;
        v_tan_norm = norm(v_tan);
        
        if v_tan_norm < 1e-8
            t_vec = [0 0 0]';
        else
            t_vec = v_tan / v_tan_norm;
        end

        [pn, pt] = Calc_Schaaf_Stress(env.rho, v_rel_mag, asin(dot(u_flow, ni)), env);
        Force_n = pn*Area_i*ni;
        Force_t = pt*Area_i*t_vec;

        F_i = Force_n + Force_t;
        
        % --- 3.5 累加力和力矩 ---
        F_hist(i,:) = F_i';
        
        % 力矩 M = r x F
        M_i = cross(r_pos(:, i), F_i);
        M_hist(i,:) = M_i';
    end
    M_hist1 = M_hist(1:n/4,:);
    M_hist2 = M_hist(n/4+1:n/2,:);
    M_hist3 = M_hist(n/2+1:3*n/4,:);
    M_hist4 = M_hist(3*n/4+1:n,:);

    M_array = [M_hist1+M_hist4; M_hist2+M_hist3];
end

function sn = sgn(x)

    if x > deg2rad(0.05)
        sn = 1;
    elseif x < -deg2rad(0.05)
        sn = -1;
    else
        sn = 0;
    end
end