% =========================================================================
% SCALES 论文绘图主程序：破损率对力矩可达域与分配精度的影响 (2x2 布局)
% =========================================================================
clear; clc; close all;

%% 1. 物理环境初始化与 a, b, c 动态提取
disp('正在加载物理仿真环境并提取解耦参数 a, b, c...');
load('para.mat'); 
sc.n = 500;
eula =[0, 0, 0]; % deg
q_err = eul2quat(deg2rad(eula), "ZYX")';

% 通过给定的 Delta 序列调用 Calculate_SCALES_plot2 提取真实力矩
Delta_init =[zeros(sc.n/2, 1); ones(sc.n/2, 1)];
M_array = Calculate_SCALES_plot2(sc, env, q_err, Delta_init);

% 提取 My 和 Mz 并计算绝对值参数
My = M_array(:, 2);
Mz = M_array(:, 3);

c_val = abs(My(1));
a_val = abs(Mz(2) - Mz(1));
b_val = abs(Mz(1)) - a_val;
m = 125;

fprintf('参数提取成功:\n a = %g\n b = %g\n c = %g\n\n', a_val, b_val, c_val);

%% 2. 蒙特卡洛打靶参数设置
num_tests = 300;  % 打靶点数（800个点能让散点图和包络线显得非常饱满美观）
target_radius = 0.06; 

% 生成统一的打靶目标池 (在半径为 0.06 的圆内均匀分布)
% rng(42); % 固定随机种子

% 1. 生成均匀分布的角度 (0 到 2*pi)
theta = 2 * pi * rand(num_tests, 1);

% 2. 生成均匀分布的半径
% 注意：使用 sqrt(rand) 是为了确保点在圆盘内是"面积均匀"的，否则中心点会过密
r = target_radius * sqrt(rand(num_tests, 1));

% 3. 转换为直角坐标系 (My, Mz)
Target_My = r .* cos(theta);
Target_Mz = r .* sin(theta);

fail_rates =[0, 0.1, 0.2, 0.4];
labels = {'(a) 0% Blades Loss', '(b) 10% Blades Loss', ...
          '(c) 20% Blades Loss', '(d) 40% Blades Loss'};

% 用于存储 4 种情况结果的 Cell 数组
True_My_all = cell(4, 1);
True_Mz_all = cell(4, 1);
Error_all   = cell(4, 1);

%% 3. 核心打靶循环
for f_idx = 1:length(fail_rates)
    p_fail = fail_rates(f_idx);
    fprintf('>>> 开始仿真子帆破损率: %d%%\n', p_fail * 100);
    
    T_My = zeros(num_tests, 1);
    T_Mz = zeros(num_tests, 1);
    Err  = zeros(num_tests, 1);
    
    tic;
    for k = 1:num_tests
        % (1) 物理层面的随机故障注入
        sail_domains = cell(4*m, 1);
        for idx = 1:4*m
            if rand() < p_fail
                sail_domains{idx} =[randi([0, 1])]; % 损坏，卡死在 0 或 1
            else
                sail_domains{idx} = [0, 1];          % 完好
            end
        end
        
        % (2) 提取宏观系流域 (根据最新物理映射规则)
        domain_u = cell(m, 1); domain_v = cell(m, 1);
        for i = 1:m
            % u = s1 - s2 (映射 1~125 & 376~500)
            dom_u1 = sail_domains{i}; dom_u2 = sail_domains{375 + i};
            u_opts =[]; for s1=dom_u1, for s2=dom_u2, u_opts(end+1)=s1-s2; end, end %#ok<*AGROW>
            domain_u{i} = unique(u_opts);
            
            % v = s2 - s1 (映射 126~250 & 251~375)
            dom_v1 = sail_domains{125 + i}; dom_v2 = sail_domains{250 + i};
            v_opts =[]; for s1=dom_v1, for s2=dom_v2, v_opts(end+1)=s2-s1; end, end 
            domain_v{i} = unique(v_opts);
        end
        
        % (3) 调用带前瞻容错的贪婪算法 (注意 My 对应常量 c 维度，Mz 对应等差维度)
        [u_calc, v_calc, ~, ~] = SCALES_FaultTolerant_Solver(a_val, b_val, c_val, Target_My(k), Target_Mz(k), domain_u, domain_v);
        
        % (4) 生成 500 位严格受限的物理控制指令
        Delta = Generate_FaultAware_Delta(u_calc, v_calc, sail_domains);
        
        % (5) 真实物理引擎解算
        M_out = Calculate_SCALES_Torque(sc, env, q_err, Delta);
        T_My(k) = M_out(2);
        T_Mz(k) = M_out(3);
        
        % 计算端到端物理误差
        Err(k) = sqrt((T_My(k) - Target_My(k))^2 + (T_Mz(k) - Target_Mz(k))^2);
    end
    toc;
    
    % 保存当前破损率的数据
    True_My_all{f_idx} = T_My;
    True_Mz_all{f_idx} = T_Mz;
    Error_all{f_idx}   = Err;
end
disp('所有蒙特卡洛仿真完成，开始绘制论文图表...');

%% 4. 论文级高质量绘图 (2x2 布局)
fig = figure('Name', 'Reachable Domain under Sub-sail Loss', ...
             'Position',[100, 100, 950, 750], 'Color', 'w');

% 使用紧凑型平铺布局 (MATLAB 2019b+)
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 寻找所有误差中的 95% 分位数，用于统一全局 Colorbar 的上限，让色彩对比更强烈
all_err_mat = cell2mat(Error_all);
max_c_limit = 0.005; % prctile(all_err_mat, 95); 

for i = 1:4
    nexttile; hold on;
    
    CX = True_My_all{i}; 
    CY = True_Mz_all{i};
    errs = Error_all{i};
    

    % [图层 1] 绘制浅灰色打靶背景池 (突出不可达区域)
    % scatter(CX, CY, 12,[0.85 0.85 0.85], 'filled');
    scatter(1000*Target_My, 1000*Target_Mz, 12,[0.85 0.85 0.85], 'filled');
    %[图层 2] 绘制实际物理落点，并根据误差热力着色
    % scatter(Target_My, Target_Mz, 18, errs, 'filled', 'MarkerEdgeColor', 'none');
    max_err = 0.01; 

    % 2. 筛选满足条件的点索引
    idx_valid = (errs <= max_err);
    
    % 3. 提取筛选后的坐标
    CX_filtered = CX(idx_valid);
    CY_filtered = CY(idx_valid);
    scatter(1000*CX_filtered, 1000*CY_filtered, 12, errs(idx_valid)*1000, 'filled', 'MarkerEdgeColor', 'none');
    % [图层 3] 绘制可达域的凸包边界 (Convex Hull) -> 视觉灵魂，展示萎缩过程
    if length(CX_filtered) >= 3
        k_hull = boundary(CX_filtered, CY_filtered, 0); % 0.7为收缩因子，能紧贴边缘
        plot(1000*CX_filtered(k_hull), 1000*CY_filtered(k_hull), 'r-', 'LineWidth', 1.5); 
    end
    
    % 坐标轴与学术排版美化
    box on; grid on;
    xlim(1000*target_radius*[-1 1]); ylim(1000*target_radius*[-1 1]);
    set(gca, 'FontSize', 10, 'FontName', 'Times New Roman', 'LineWidth', 0.8, 'Layer', 'top');
    
    % 设置颜色映射和统一色轴限制
    colormap("jet");
    caxis([0, max_c_limit*1000]); % 如果是极新版 MATLAB 可替换为 clim([0, max_c_limit])
    
    % 添加 X/Y 轴标签 (仅在外部边缘显示，保持图面干净)
    if i == 3 || i == 4
        xlabel('$\it{M_y}$ $\rm (mN \cdot m$)', 'Interpreter', 'latex', 'FontSize', 10);
    end
    if i == 1 || i == 3
        ylabel('$\it{M_z}$ $\rm (mN \cdot m$)', 'Interpreter', 'latex', 'FontSize', 10);
    end
    
    title(labels{i}, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'normal');
end

% 增加全局 Colorbar
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Allocation Error (mN\cdotm)';
cb.Label.Interpreter = 'tex';
cb.Label.FontName = 'Times New Roman';
cb.Label.FontSize = 10;

% 增加图例解释说明
lgd = legend('Target Torque', 'Achieved Torque', 'Reachable Boundary', ...
             'Orientation', 'horizontal', 'NumColumns', 3);
lgd.Layout.Tile = 'north';
lgd.FontName = 'Times New Roman';
lgd.FontSize = 10;

disp('绘图完成！建议在 Figure 窗口点击 "文件 -> 导出设置" 以 600dpi 导出为 TIFF 或 EMF 格式插入 Word。');


function Delta = Generate_FaultAware_Delta(u, v, sail_domains)
    % 严格遵循坏点约束，生成能够输出目标 u,v 的真实子帆控制指令 Delta
    m = 125;
    Delta = zeros(4*m, 1);
    
    for i = 1:m
        % ---- 处理 u 序列 (映射到 i 和 375+i) ----
        p1 = i; p2 = 375 + i;
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
        
        % ---- 处理 v 序列 (映射到 125+i 和 250+i) ----
        % 注意: v 的逆向规则 v=1 对应 0 1, v=0 对应 1 1, v=-1 对应 1 0
        p1 = 125 + i; p2 = 250 + i;
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

function dom = get_sub_sail_domain(prob_fail)
    % 生成单个子帆面合法的状态集合
    if rand() < prob_fail
        dom = [randi([0, 1])]; % 发生卡死，固定在 0 或 1
    else
        dom = [0, 1];          % 健康，可取 0 或 1
    end
end