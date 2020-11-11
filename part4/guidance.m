function [ chi_d ] = guidance(x_t, y_t, x_ref, y_ref, x, y,L, u,v)
    pi_p = atan2(y_t-y_ref,x_t-x_ref);
    %xhi_inf = pi/2;
    [y_e] = crosstrackWpt(x_t, y_t, x_ref, y_ref, x, y);
    delta = 10*L;
    K_p = 1/delta; %2*L
    y_in
    chi_d = pi_p - atan(K_p*y_e); %- asin(v/(sqrt(u^2+v^2)));
    chi_d = wrapTo2Pi(chi_d);
    %chi_d = ssa(chi_d);
end

