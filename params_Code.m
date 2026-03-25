% params

%% 1. 系统与环境初始化
% close all
% clear
% --- SCALES 参数 ---
dens_sail = 3.11e-3;
dens_frame = 0.0961;

sc.s  = 0.3;

sc.lb = 3 ;

l = 5 + sc.s;
sc.wb = (l - sc.s)/sc.n*2; % 适当调整尺寸以匹配 n
w = sc.lb*2;
sc.m1 = 20; % dens_kern*(sc.s)^3;
sc.m2 = 3e-3;
sc.m3 = dens_frame*(l + w) + 2*dens_frame*(l + w);
sc.m4 = sc.lb*sc.wb*dens_sail;

sc.l = l;
sc.w = w;

J_b = Calculate_SCALES_Inertia(sc);
sc.J = J_b;

% --- 环境参数 (350km LEO) ---
env.h = 280e3; env.rho = 5.24e-11;
env.f107 = 150;
env.Tw = 300; env.Tinf = 1000; env.M = 0.016;
env.sigma_n = 0.87; env.sigma_t = 0.87;

% save("para.mat")