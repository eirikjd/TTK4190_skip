function [ chi_d ] = guidance(x_t, y_t, x_ref, y_ref, x, y,L)
    pi_p = atan2(y_t-y_ref,x_t-x_ref);
    %xhi_inf = pi/2;
    [x_p, y_p, y_e] = crosstrack(x_t, y_t, x_ref, y_ref, x, y);
    K_p = 200; %2*L
    chi_d = pi_p - arctan(K_p*y_e);
    chi_d = ssa(chi_d);
end

