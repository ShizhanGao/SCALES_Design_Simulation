function[ans_u, ans_v, final_X, final_Y] = SCALES_FaultTolerant_Solver(a, b, c, target_X, target_Y, domain_u, domain_v)
% SCALES_FAULTTOLERANT_SOLVER 带有底层子帆容错机制的贪婪微调求解器
%
% 输入参数:
%   a, b, c    : 系统的线性解耦参数
%   target_X   : 常量维度目标 (Z维度)
%   target_Y   : 等差维度目标 (W维度)
%   domain_u   : 1×m 的 cell 数组，每个元素是该点 u_i 允许的取值集合 (如[-1, 0])
%   domain_v   : 1×m 的 cell 数组，每个元素是该点 v_i 允许的取值集合 (如 [0, 1])
%
% 输出参数:
%   ans_u, ans_v : 最终分配的系数 (-1, 0, 1)
%   final_X, Y   : 算法最终合成的实际坐标

    m = length(domain_u);
    A_vals = b + (1:m)' * a;
    Z_req = round(target_X / c);
    
    %% 第一阶段：预处理与局部动态菜单生成 (Pre-processing)
    % 我们将 W 映射到 1~5 的索引 (W + 3) 以极速替代 Hash Map
    valid_W = cell(m, 1);
    valid_Z = cell(m, 5);
    min_Z = nan(m, 5); 
    max_Z = nan(m, 5);
    
    global_min_Z = zeros(m, 1); % 每个位置能提供的绝对最小 Z
    global_max_Z = zeros(m, 1); % 每个位置能提供的绝对最大 Z
    
    for i = 1:m
        u_opts = domain_u{i};
        v_opts = domain_v{i};
        
        % 穷举当前位置所有合法的 (W, Z) 组合
        w_list = []; z_list =[];
        for u = u_opts
            for v = v_opts
                w_list(end+1) = u + v; %#ok<AGROW>
                z_list(end+1) = u - v; %#ok<AGROW>
            end
        end
        
        unique_w = unique(w_list);
        valid_W{i} = unique_w;
        
        curr_min = inf; curr_max = -inf;
        for w = unique_w
            idx = w + 3;
            zs = unique(z_list(w_list == w));
            valid_Z{i, idx} = zs;
            min_Z(i, idx) = min(zs);
            max_Z(i, idx) = max(zs);
            
            curr_min = min(curr_min, min(zs));
            curr_max = max(curr_max, max(zs));
        end
        global_min_Z(i) = curr_min;
        global_max_Z(i) = curr_max;
    end
    
    % 计算"前瞻包络"（未来 1~i-1 还能提供的极限容量）
    fut_min_Z =[0; cumsum(global_min_Z(1:end-1))]; 
    fut_max_Z =[0; cumsum(global_max_Z(1:end-1))];

    %% 第二阶段：带前瞻包络保护的贪婪搜索 (Greedy with Look-ahead)
    W_sel = zeros(m, 1);
    curr_min_Z = 0; curr_max_Z = 0;
    rem_Y = target_Y;
    
    for i = m:-1:1
        w_opts = valid_W{i};
        % 按逼近 Y 的效果对合法 W 进行排序
        errs = abs(rem_Y - w_opts * A_vals(i));
        [~, sort_idx] = sort(errs);
        w_opts_sorted = w_opts(sort_idx);
        
        best_w = w_opts_sorted(1); % 默认兜底
        for w_test = w_opts_sorted
            z_min_test = min_Z(i, w_test + 3);
            z_max_test = max_Z(i, w_test + 3);
            
            % 核心包络监控：选了这个 w，加上未来的极限容量，还能包住 Z_req 吗？
            env_min = curr_min_Z + z_min_test + fut_min_Z(i);
            env_max = curr_max_Z + z_max_test + fut_max_Z(i);
            
            if Z_req >= env_min && Z_req <= env_max
                best_w = w_test; break; % 完美包络，就决定是它了
            end
        end
        
        W_sel(i) = best_w;
        curr_min_Z = curr_min_Z + min_Z(i, best_w + 3);
        curr_max_Z = curr_max_Z + max_Z(i, best_w + 3);
        rem_Y = rem_Y - best_w * A_vals(i);
    end
    
    %% 第三阶段：奇偶对齐修复 (Parity Alignment)
    % 确保 sum(W) 与 Z_req 奇偶性一致
    if mod(sum(W_sel), 2) ~= mod(Z_req, 2)
        best_k = 0; best_w_new = 0; min_err = inf;
        best_dmin = 0; best_dmax = 0;
        
        for k = 1:m
            for w_new = valid_W{k}
                if abs(w_new - W_sel(k)) == 1 % 寻找改变量为1的替换方案
                    d_min = min_Z(k, w_new+3) - min_Z(k, W_sel(k)+3);
                    d_max = max_Z(k, w_new+3) - max_Z(k, W_sel(k)+3);
                    % 校验替换后包络是否依然安全
                    if Z_req >= curr_min_Z + d_min && Z_req <= curr_max_Z + d_max
                        err = abs(rem_Y - (w_new - W_sel(k)) * A_vals(k));
                        if err < min_err
                            min_err = err; best_k = k; best_w_new = w_new;
                            best_dmin = d_min; best_dmax = d_max;
                        end
                    end
                end
            end
        end
        if best_k > 0
            rem_Y = rem_Y - (best_w_new - W_sel(best_k)) * A_vals(best_k);
            W_sel(best_k) = best_w_new;
            curr_min_Z = curr_min_Z + best_dmin;
            curr_max_Z = curr_max_Z + best_dmax;
        end
    end

    %% 第四阶段：菜单约束下的连续平移微调 (Domain-Aware Fine-Tuning)
    if a ~= 0
        for iter = 1:20
            delta_steps = round(rem_Y / a);
            if delta_steps == 0; break; end
            applied = false;
            
            direction = sign(delta_steps);
            steps = abs(delta_steps);
            
            for gap = min(steps, m-1):-1:1
                for i = 1 : m - gap
                    j = i + gap;
                    w_i_new = W_sel(i) - direction;
                    w_j_new = W_sel(j) + direction;
                    
                    % 严苛的合法性校验：减前加后是否都在彼此的菜单里？
                    if ismember(w_i_new, valid_W{i}) && ismember(w_j_new, valid_W{j})
                        d_min = min_Z(i, w_i_new+3) - min_Z(i, W_sel(i)+3) + ...
                                min_Z(j, w_j_new+3) - min_Z(j, W_sel(j)+3);
                        d_max = max_Z(i, w_i_new+3) - max_Z(i, W_sel(i)+3) + ...
                                max_Z(j, w_j_new+3) - max_Z(j, W_sel(j)+3);
                        
                        % 校验包络
                        if Z_req >= curr_min_Z + d_min && Z_req <= curr_max_Z + d_max
                            W_sel(i) = w_i_new; W_sel(j) = w_j_new;
                            curr_min_Z = curr_min_Z + d_min;
                            curr_max_Z = curr_max_Z + d_max;
                            rem_Y = rem_Y - direction * (A_vals(j) - A_vals(i));
                            applied = true; break;
                        end
                    end
                end
                if applied; break; end
            end
            if ~applied; break; end
        end
    end

    %% 第五阶段：菜单内的精确对消吸收 (Exact Z Matching)
    Z_sel = zeros(m, 1);
    for i = 1:m
        Z_sel(i) = min_Z(i, W_sel(i) + 3); % 初始取下界
    end
    diff_Z = Z_req - sum(Z_sel);
    
    % 因为奇偶已对齐，diff_Z 必定是偶数，且只能增加
    while diff_Z > 0
        applied = false;
        for i = 1:m
            zs = valid_Z{i, W_sel(i) + 3};
            % 因为物理对称性，合法的 Z 必然按步长 2 递增
            if Z_sel(i) + 2 <= max(zs) && ismember(Z_sel(i) + 2, zs)
                Z_sel(i) = Z_sel(i) + 2;
                diff_Z = diff_Z - 2;
                applied = true;
                if diff_Z == 0; break; end
            end
        end
        if ~applied; break; end % 安全逃逸
    end

    %% 最终反算真实系数
    ans_u = (W_sel + Z_sel) / 2;
    ans_v = (W_sel - Z_sel) / 2;
    
    final_X = c * sum(ans_u - ans_v);
    final_Y = sum((ans_u + ans_v) .* A_vals);
end