function [pn, pt] = Calc_Schaaf_Stress(rho, V_rel, alpha_rad, params)
% CALC_SCHAAF_STRESS 计算稀薄气体法向压强和切向应力
%
% 输入:
%   rho       : 大气密度 (kg/m^3)
%   V_rel     : 相对速度大小 (m/s)
%   alpha_rad : 攻角 (角度制, 0=平行气流, pi/2=垂直气流)
%   params    : 包含物理参数的结构体 (.Tw, .Tinf, .M, .sigma)
%
% 输出:
%   pn        : 法向压强 (Pascal, N/m^2)
%   pt        : 切向/剪切应力 (Pascal, N/m^2)

    %% 1. 提取参数与基础计算
    Tw = params.Tw;       % 壁面温度 (K)
    Tinf = params.Tinf;   % 来流大气温度 (K)
    M = params.M;         % 气体摩尔质量 (kg/mol)
    sn = params.sigma_n;  % 法向动量适应系数
    st = params.sigma_t;  % 切向动量适应系数
    
    R_gas = 8.314;        % 通用气体常数
    
    % 气体常数 (J/kg*K)
    R_spec = R_gas / M;   
    
    % 最概然热运动速度 (Most Probable Thermal Speed)
    V_th = sqrt(2 * R_spec * Tinf);
    
    % 分子速度比 (Speed Ratio)
    S = V_rel / V_th;
    
    % 几何参数 s (分子速度比在法向的投影)
    % alpha 为流与板夹角: sin(alpha) 对应法向分量
    s = S * sin(alpha_rad); 
    
    %% 2. Schaaf & Chambre 公式实现
    sqrt_pi = sqrt(pi);
    Tr_ratio = sqrt(Tw / Tinf); % 温度比项
    
    % 预计算指数和误差函数项
    exp_s2 = exp(-s.^2);
    erf_s  = erf(s);
    
    % --- 计算法向压强 pn ---
    % 公式结构: 动压/S^2 * [ (镜面反射项) + (漫反射项) ]
    term_n1 = (2 - sn)/sqrt_pi * s + sn/2 * Tr_ratio;
    term_n2 = (2 - sn)*(s.^2 + 0.5) + sn/2 * Tr_ratio * sqrt_pi * s;
    
    % 系数 Cn' (注意这不是无量纲系数，包含了分母处理)
    % 这里的 pn 直接是应力 (Pa)
    coeff_common = (0.5 * rho * V_rel^2) / (S^2);
    
    pn = coeff_common * ( term_n1 * exp_s2 + term_n2 * (1 + erf_s) );
    
    % --- 计算切向应力 pt ---
    % 注意切向投影: cos(alpha)
    term_t = exp_s2 + sqrt_pi * s * (1 + erf_s);
    
    % 切向应力公式
    pt = coeff_common * (st * S * cos(alpha_rad) / sqrt_pi) * term_t;

end