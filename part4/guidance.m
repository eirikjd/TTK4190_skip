function [ chi_d ] = guidance(x_t, y_t, x_ref, y_ref, x, y)
    pi_p = atan2(y_t-y_ref,x_t-x_ref);
    %xhi_inf = pi/2;
    [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y);
    chi_d = pi_p - arctan(K_p_wp*y_e);
end

