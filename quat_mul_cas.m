function qout = quat_mul_cas(a, b)
    % Hamilton product for 4x1 quaternions (a and b are 4x1)
    % a = [a0; av], b = [b0; bv]
    a0 = a(1); av = a(2:4);
    b0 = b(1); bv = b(2:4);
    q0 = a0*b0 - av' * bv;
    qv = a0 * bv + b0 * av + cross(av, bv);
    qout = [q0; qv];
end