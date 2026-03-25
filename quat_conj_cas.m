function qc = quat_conj_cas(q)
    % q: 4x1 (SX or numeric), return 4x1 conj
    qc = [q(1); -q(2); -q(3); -q(4)];
end