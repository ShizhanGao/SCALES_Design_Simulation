function J_b = Calculate_SCALES_Inertia(para)
% CALCULATE_SCALES_INERTIA 计算分布式大气帆(SCALES)的转动惯量矩阵
%
% 输入参数 para 是一个结构体，必须包含以下字段:
%   para.m1  - 卫星核心(Kernel)质量 [kg]
%   para.m2  - 单个控制单元(Control Unit)质量 [kg]
%   para.m3  - 框架(Frame)总质量 [kg]
%   para.m4  - 单个子帆(Blade)质量 [kg]
%   para.s   - 卫星核心立方体边长 [m]
%   para.n   - 子帆总数量 (必须是4的倍数)
%   para.lb  - 单个子帆长度 (对应框架高度的一半) [m]
%   para.wb  - 单个子帆宽度 [m]
%
% 输出参数:
%   J_b      - 3x3 对角转动惯量矩阵 diag([Jx, Jy, Jz]) [kg*m^2]
%
% 坐标系定义:
%   x轴: 垂直于帆面
%   y轴: 沿帆面长度方向 (水平)
%   z轴: 沿帆面宽度方向 (竖直)

    %% 1. 提取参数与几何计算
    m1 = para.m1;
    m2 = para.m2;
    m3 = para.m3;
    m4 = para.m4;
    s  = para.s;
    n  = para.n;
    lb = para.lb;
    wb = para.wb;

    % 框架几何尺寸近似
    % 框架总长 l (由n/2个子帆宽度加上核心边长组成)
    l = (n / 2) * wb + s; 
    % 框架总高 w (由上下两个子帆长度组成)
    % w = 2 * lb; % 代码中直接使用 lb 计算，无需显式定义 w

    %% 2. 计算累加项 A_sum
    % A_sum 代表所有分布在y轴上的组件(子帆和控制单元)相对于原点的距离平方和
    % 逻辑：左右对称，共有 n/2 列；考虑核心偏移 s
    
    % 每一侧的列数 (假设左右对称，每列上下各一个，故列数为 n/2，单侧 n/4)
    num_cols_per_side = n / 4;
    
    % 使用向量化计算代替循环
    k = 1:num_cols_per_side; 
    
    % 第 k 列子帆中心距离原点的 y 距离公式推导:
    % dist = s/2 + (k-1)*wb + wb/2 = 0.5 * (s + (2*k - 1)*wb)
    % 这里直接计算累加项公式: sum( (s + (2k-1)wb)^2 )
    % 注意：物理上总方差 = 4 * sum(y_k^2) = 4 * 0.25 * sum(...) = sum(...)
    term_k = s + (2.*k - 1) * wb;
    A_sum = sum(term_k.^2);

    %% 3. 分项转动惯量计算
    
    % --- 核心 (Kernel) ---
    J_kernel = (1/6) * m1 * s^2; % Jx=Jy=Jz for cube

    % --- 框架 (Frame) ---
    % 绕 Y 轴
    J_frame_y = (2 * m3 * lb^2 * (l + lb)) / (3 * (l + 2*lb));
    % 绕 Z 轴
    J_frame_z = (m3 * l^2 * (l + 4*lb)) / (12 * (l + 2*lb));

    % --- 组合计算 Jy (绕水平轴) ---
    % 组成: 核心 + n个控制单元(平移) + n个子帆(自转+平移) + 框架
    % 子帆绕自身水平轴: 1/12*m4*lb^2, 平移 lb/2 -> 1/3*m4*lb^2
    J_y = J_kernel + ...
          n * m2 * lb^2 + ...
          (1/3) * n * m4 * lb^2 + ...
          J_frame_y;

    % --- 组合计算 Jz (绕竖直轴) ---
    % 组成: 核心 + (控制单元+子帆)平移 + 子帆自转 + 框架
    % 平移项利用 A_sum
    J_z = J_kernel + ...
          (m2 + m4) * A_sum + ...
          (n * m4 * wb^2) / 12 + ...
          J_frame_z;

    % --- 组合计算 Jx (垂直帆面轴) ---
    % 利用垂直轴定理扩展: Jx = Jy + Jz - Z轴方向重叠的核心惯量
    % 因为 J_kernel 在 Jy 和 Jz 中都加了一次，对于立体核心，
    % Jx_kernel = Jy_kernel = Jz_kernel
    % 帆面部分满足 Jx_sail = Jy_sail + Jz_sail
    % 所以 Jx = (Jy - J_kernel) + (Jz - J_kernel) + J_kernel
    J_x = J_y + J_z - J_kernel;

    %% 4. 组装输出矩阵
    J_b = diag([J_x, J_y, J_z]);

end