%==========================================================================
% Collagen part of strain energy, derivative w.r.t prestretch x
function y=dWcdx(x,kc)
% --------
% by Hamid, 2018
%==========================================================================
    y=kc(2)*(x*(x^2-1))*exp(kc(3)*(x^2-1)^2);
end

