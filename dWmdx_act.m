%==========================================================================
% SMC active part of strain energy, derivative w.r.t prestretch
function y=dWmdx_act(x)
% --------
% by Hamid, 2018
%==========================================================================
    global rho_w S lmax lmin

    lm=lmax;
    l0=lmin;

    y=S/rho_w*(1-(lm-x)^2/(lm-l0)^2);
end

