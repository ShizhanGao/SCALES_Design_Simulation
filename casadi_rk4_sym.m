function x_next = casadi_rk4_sym(x, u, h, I, Iinv)
    k1 = dynamics_cas(x, u, I, Iinv);
    k2 = dynamics_cas(x + 0.5*h*k1, u, I, Iinv);
    k3 = dynamics_cas(x + 0.5*h*k2, u, I, Iinv);
    k4 = dynamics_cas(x + h*k3, u, I, Iinv);
    x_next = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end