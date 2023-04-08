%==========================================================================
% SMC passive part of strain energy, derivative w.r.t prestretch
function y=dWmdx(x,kc)
% --------
% by Hamid, 2018
%==========================================================================
    y=kc(4)*(x*(x^2-1))*exp(kc(5)*(x^2-1)^2);
end

