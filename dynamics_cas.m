function xdot = dynamics_cas(x, u, I, Iinv)
    % x: [q0 q1 q2 q3 wx wy wz] as CasADi SX column
    q = x(1:4);
    w = x(5:7);
    we = w - [0; 1.1455e-3; 0];
    % Omega matrix (symbolic)
    Omega = [    0    -we(1) -we(2) -we(3);
              we(1)     0    we(3) -we(2);
              we(2) -w(3)     0   we(1);
              we(3)  we(2) -we(1)     0 ];
    qdot = 0.5 * (Omega * q);
   
    %% 
    C3 = [2*(q(2)*q(4) - q(1)*q(3));
          2*(q(3)*q(4) + q(1)*q(2));
          1 - 2*(q(2)^2 + q(3)^2)];
    Mgg = 3*1.1455e-3^2*[(I(3) - I(2))*C3(2)*C3(1);
                         (I(1) - I(3))*C3(2)*C3(1);
                         (I(2) - I(1))*C3(2)*C3(1)];
    %%
    wdot = Iinv * (u - cross(w, I*w) + Mgg);

    xdot = [qdot; wdot];
end