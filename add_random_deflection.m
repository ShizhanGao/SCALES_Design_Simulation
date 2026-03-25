function v_deflected = add_random_deflection(v_original, angle_deg)
    % ADD_RANDOM_DEFLECTION 对矢量施加随机的定值偏转
    % 输入:
    %   v_original: 原始矢量 (3x1 或 1x3)
    %   angle_deg:  偏转角度 (度)，例如 3
    % 输出:
    %   v_deflected: 偏转后的矢量 (保持模长不变)

    % 1. 预处理
    v_col = v_original(:); % 确保是列向量
    mag = norm(v_col);     % 获取模长
    
    % if mag < 1e-12
    %     error('零向量无法偏转方向');
    % end
    
    n_vec = v_col / mag;   % 归一化主轴向量
    
    % 2. 构建局部正交基 (Local Basis)
    % null(n_vec') 会返回两个互相垂直且都垂直于 n_vec 的单位向量
    % 构成了零空间基，正好可以用作 b1 和 b2
    basis = null(n_vec'); 
    b1 = basis(:, 1);
    b2 = basis(:, 2);
    
    % 3. 生成随机方位角 phi
    phi = 2 * pi * rand();
    theta = deg2rad(angle_deg);
    
    % 4. 合成偏转后的单位矢量
    % 公式: v' = n*cos(theta) + b1*sin(theta)cos(phi) + b2*sin(theta)sin(phi)
    v_deflected_unit = n_vec * cos(theta) + ...
                       b1 * (sin(theta) * cos(phi)) + ...
                       b2 * (sin(theta) * sin(phi));
    
    % 5. 恢复模长
    v_deflected = v_deflected_unit * mag;
    
    % (可选) 如果输入是行向量，输出也转为行向量
    if isrow(v_original)
        v_deflected = v_deflected';
    end
end
