function solver = build_ancillary_mpc(Ir, dt, N)
    
    import casadi.*
    nx = 7; nu = 2;
    
    % --- 1. 决策变量 ---
    u_anc = SX.sym('u_anc', nu*N); % [u_0; u_1; ...; u_{N-1}]
    
    % --- 2. 参数 p 的定义 (关键修正) ---
    % p 由两部分组成:
    % Part 1: 初始状态 x0 (7维) % 真实轨迹0时刻状态
    % Part 2: 参考轨迹序列 ref_seq. 
    %         我们需要 0~N-1 时刻的参考值来算过程Cost，
    %         需要 N 时刻的参考值来算终端Cost/Constraint。
    %         总共 N+1 个块。每个块包含 [z_ref(7); u_ref(2)]。
    
    x_init_sym = SX.sym('x_init', nx);
    
    len_blk = nx + nu; % 每个参考块的长度 = 9
    ref_seq_sym = SX.sym('ref_seq', len_blk * (N+1)); 
    
    % 将它们拼接成唯一的参数向量 p
    p_vec = [x_init_sym; ref_seq_sym];
    
    % --- 3. 权重定义 ---
    Q_q = 100;    % 过程姿态权重
    Q_w = 10;     % 过程角速度权重
    R_u = 0.1;    % 控制权重
    P_q = 1000;   % 终端姿态权重 (Terminal Cost, lambda*Vf)
    P_w = 100;    % 终端角速度权重
    
    % 转换转动惯量
    I_c = SX(Ir); 
    Iinv_c = SX(inv(Ir)); 
    I3 = SX.eye(3);

    cost = 0;
    xk = x_init_sym; % 状态变量初始化为 p 中的初始状态
    
    % ==========================================================
    % Loop: 预测 0 到 N-1 步 (Stage Cost & Integration)
    % ==========================================================
    for k = 0:N-1
        % --- A. 从参数 p 中提取第 k 步的参考值 ---
        % ref_seq_sym 的结构: [Block_0; Block_1; ...; Block_N]
        % Matlab索引从1开始
        idx_start = k * len_blk + 1; 
        
        % 提取 z_ref_k (7维) 和 u_ref_k (2维)
        z_ref_k = ref_seq_sym(idx_start : idx_start + nx - 1);
        u_ref_k = ref_seq_sym(idx_start + nx : idx_start + len_blk - 1);
        
        q_ref = z_ref_k(1:4);
        w_ref = z_ref_k(5:7);
        
        % --- B. 从决策变量中提取第 k 步控制 ---
        uk = u_anc(k*nu + 1 : (k+1)*nu);
        
        % --- C. 计算 Stage Cost (l(x,u)) ---
        % 1. 四元数误差
        qe = quat_mul_cas(quat_conj_cas(q_ref), xk(1:4));
        e0 = qe(1); e_vec = qe(2:4);
        
        % 2. 角速度误差 (带旋转矩阵修正)
        skew_e = [0, -e_vec(3), e_vec(2); e_vec(3), 0, -e_vec(1); -e_vec(2), e_vec(1), 0];
        R_mat = (e0^2 - dot(e_vec, e_vec))*I3 + 2*(e_vec*e_vec') - 2*e0*skew_e;
        we = xk(5:7) - mtimes(R_mat, w_ref);
        
        % 3. 控制误差
        ue = uk - u_ref_k;
        
        step_cost = Q_q * (1 - abs(e0))^2 + ...
                    Q_w * (we' * we) + ...
                    R_u * (ue' * ue);
        cost = cost + step_cost;
        
        % --- D. 状态积分 (rk4) ---
        ufull = [SX(0); uk]; % 欠驱动
        xk = casadi_rk4_sym(xk, ufull, dt, I_c, Iinv_c);
        
        % 归一化
        xk(1:4) = xk(1:4) / sqrt(sumsqr(xk(1:4)));
    end
    
    % ==========================================================
    % Terminal Step: 第 N 步 (Terminal Cost & Constraint)
    % ==========================================================
    % 此时 xk 已经是 x_N (预测末端状态)
    
    % 1. 从参数 p 中提取第 N 步的参考值
    idx_final = N * len_blk + 1;
    z_ref_N = ref_seq_sym(idx_final : idx_final + nx - 1);
    % 注意: u_ref_N 通常设为0或不需要，这里只需要 z_ref_N
    
    q_ref_N = z_ref_N(1:4);
    w_ref_N = z_ref_N(5:7);
    
    % 2. 计算末端误差
    qe_N = quat_mul_cas(quat_conj_cas(q_ref_N), xk(1:4));
    e0_N = qe_N(1); e_vec_N = qe_N(2:4);
    
    skew_e_N = [0, -e_vec_N(3), e_vec_N(2); e_vec_N(3), 0, -e_vec_N(1); -e_vec_N(2), e_vec_N(1), 0];
    R_mat_N = (e0_N^2 - dot(e_vec_N, e_vec_N))*I3 + 2*(e_vec_N*e_vec_N') - 2*e0_N*skew_e_N;
    we_N = xk(5:7) - mtimes(R_mat_N, w_ref_N);
    
    % --- E. 添加 Terminal Cost (lambda * Vf) ---
    % 这是一个非常重要的修正，保证稳定性
    term_cost = P_q * (1 - abs(e0_N))^2 + P_w * (we_N' * we_N);
    cost = cost + term_cost;
    
    % --- F. 定义 Terminal Constraint (g) ---
    % g = [qe_scalar; qe_vector; w_error]
    % 对应 bounds: [0.995; -eps; -eps...; -eps] <= g <= [1; eps; ... eps]
    g_term = [e0_N; e_vec_N; we_N];
    
    % --- 构建 NLP ---
    nlp = struct('x', u_anc, 'f', cost, 'g', g_term, 'p', p_vec);
    
    opts = struct;
    opts.ipopt.print_level = 0; 
    opts.print_time = 0;
    opts.ipopt.max_iter = 500;
    opts.ipopt.tol = 1e-8;
    
    solver = nlpsol('ancillary_solver', 'ipopt', nlp, opts);
end