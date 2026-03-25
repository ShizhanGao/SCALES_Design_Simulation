function x1 = rk4_step_numeric(x, u, h, I, Iinv)
    k1 = dynamics_numeric(x, u, I, Iinv);
    k2 = dynamics_numeric(x + 0.5*h*k1, u, I, Iinv);
    k3 = dynamics_numeric(x + 0.5*h*k2, u, I, Iinv);
    k4 = dynamics_numeric(x + h*k3, u, I, Iinv);
    x1 = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end