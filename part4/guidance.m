function [ chi_d ] = guidance(x_t, y_t, x_ref, y_ref, x, y, L)
    pi_p = atan2(y_t-y_ref,x_t-x_ref);
    %xhi_inf = pi/2;
    K_p_wp = 200;%1/(5 * L);
    [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y);
    chi_d = pi_p - atan(K_p_wp*y_e);
    chi_d = ssa(chi_d);
end

