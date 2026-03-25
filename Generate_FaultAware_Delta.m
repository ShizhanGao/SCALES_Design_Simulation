function Delta = Generate_FaultAware_Delta(u, v, sail_domains)
    % 严格遵循坏点约束，生成能够输出目标 u,v 的真实子帆控制指令 Delta
    m = length(u);
    Delta = zeros(4*m, 1);
    
    for i = 1:m
        % ---- 处理 u 序列 (映射到 i 和 3*m+i) ----
        p1 = i; p2 = 3*m + i;
        dom1 = sail_domains{p1}; dom2 = sail_domains{p2};
        if u(i) == 1 && ismember(1, dom1) && ismember(0, dom2)
            Delta(p1) = 1; Delta(p2) = 0;
        elseif u(i) == -1 && ismember(0, dom1) && ismember(1, dom2)
            Delta(p1) = 0; Delta(p2) = 1;
        elseif u(i) == 0 
            if ismember(1, dom1) && ismember(1, dom2)
                Delta(p1) = 1; Delta(p2) = 1; % 优先选用默认的 1 1 状态
            elseif ismember(0, dom1) && ismember(0, dom2)
                Delta(p1) = 0; Delta(p2) = 0; % 坏点迫使使用 0 0 状态
            else
                error('物理映射异常：u=0 无解');
            end
        else
            error('物理映射异常：u 不合规');
        end
        
        % ---- 处理 v 序列 (映射到 m+i 和 2*m+i) ----
        % 注意: v 的逆向规则 v=1 对应 0 1, v=0 对应 1 1, v=-1 对应 1 0
        p1 = m + i; p2 = 2*m + i;
        dom1 = sail_domains{p1}; dom2 = sail_domains{p2};
        if v(i) == 1 && ismember(0, dom1) && ismember(1, dom2)
            Delta(p1) = 0; Delta(p2) = 1;
        elseif v(i) == -1 && ismember(1, dom1) && ismember(0, dom2)
            Delta(p1) = 1; Delta(p2) = 0;
        elseif v(i) == 0 
            if ismember(1, dom1) && ismember(1, dom2)
                Delta(p1) = 1; Delta(p2) = 1;
            elseif ismember(0, dom1) && ismember(0, dom2)
                Delta(p1) = 0; Delta(p2) = 0;
            else
                error('物理映射异常：v=0 无解');
            end
        else
            error('物理映射异常：v 不合规');
        end
    end
end