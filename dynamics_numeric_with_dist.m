function dy = dynamics_numeric_with_dist(t, y, u, It, dist_torque)
    q = y(1:4);
    w = y(5:7);
    
    % Kinematics
    G = [-q(2) -q(3) -q(4);
          q(1) -q(4)  q(3);
          q(4)  q(1) -q(2);
         -q(3)  q(2)  q(1)];
    we = w - [0; 1.1455e-3; 0];
    dq = 0.5 * G * we;
    
    % Dynamics (Euler) with Disturbance
    % It * dw + w x (It*w) = u + dist
    C3 = [2*(q(2)*q(4) - q(1)*q(3));
          2*(q(3)*q(4) + q(1)*q(2));
          1 - 2*(q(2)^2 + q(3)^2)];
    Mgg = 3*1.1455e-3^2*[(It(3) - It(2))*C3(2)*C3(1);
                         (It(1) - It(3))*C3(2)*C3(1);
                         (It(2) - It(1))*C3(2)*C3(1)];

    u_full = u + dist_torque + Mgg; 
    
    dw = It \ (u_full - cross(w, It*w));
    
    dy = [dq; dw];
end