function [z_ref, u_ref] = interpolate_reference(t_query, t_ref_arr, Xhist, Uhist)
    % t_query: 当前需要查询的时间点
    % t_ref_arr: 参考轨迹的时间轴 (0:dtr:TotalTime)
    % Xhist: 7x(Nr+1)
    % Uhist: 2xNr
    
    nx = 7; nu = 2;
    
    % 1. 处理时间越界 (如果预测时间超过了参考轨迹总长，保持在最后状态)
    if t_query >= t_ref_arr(end)
        z_ref = Xhist(:, end);
        u_ref = zeros(nu, 1); % 假设到达终点后不再施加前馈控制
        return;
    end
    
    % 2. 查找索引
    % 找到 idx 使得 t_ref_arr(idx) <= t_query < t_ref_arr(idx+1)
    idx = find(t_ref_arr <= t_query, 1, 'last');
    if isempty(idx), idx = 1; end
    if idx >= length(t_ref_arr), idx = length(t_ref_arr) - 1; end
    
    % 3. 计算线性插值系数 alpha
    t0 = t_ref_arr(idx);
    t1 = t_ref_arr(idx+1);
    alpha = (t_query - t0) / (t1 - t0);
    
    % 4. 状态插值
    q0 = Xhist(1:4, idx); q1 = Xhist(1:4, idx+1);
    w0 = Xhist(5:7, idx); w1 = Xhist(5:7, idx+1);
    
    % 四元数插值 (归一化线性插值，近似 Slerp)
    q_interp = (1 - alpha) * q0 + alpha * q1;
    q_interp = q_interp / norm(q_interp);
    
    % 角速度插值
    w_interp = (1 - alpha) * w0 + alpha * w1;
    z_ref = [q_interp; w_interp];
    
    % %% 5. 控制量插值
    % u0 = Uhist(:, idx); 
    % % 注意 Uhist 长度比 Xhist 少 1，如果 idx 是最后一个点需特殊处理
    % if idx + 1 > size(Uhist, 2)
    %     u1 = zeros(nu, 1);
    % else
    %     u1 = Uhist(:, idx+1);
    % end
    % u_ref = (1 - alpha) * u0 + alpha * u1;

    %% 5. 控制量提取 (修改为零阶保持 Zero-Order Hold)
    % 逻辑：当 t 位于 [t_ref(idx), t_ref(idx+1)) 区间时，
    % 控制量应恒定为 Uhist(:, idx)，不需要 alpha 参与插值。
    
    if idx > size(Uhist, 2)
        % 边界情况：如果查询时间点正好是轨迹结束时刻或超出范围
        % 此时没有定义的控制量，通常设为 0
        u_ref = zeros(nu, 1);
    else
        % 直接取当前区间的控制量
        u_ref = Uhist(:, idx);
    end
end