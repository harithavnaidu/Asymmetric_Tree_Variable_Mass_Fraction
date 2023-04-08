%==========================================================================
% Get flow from Hfn, frequency domain
function [Qn] = Qrecovery(NumModes,omegan,L,x,Hfn,cn,Z0,Gamma)
% ---------
% by Vasilina, 2018
%==========================================================================   
    Qn = zeros(NumModes,1);
    for k=2:NumModes
        Qn(k) = Hfn(k) *exp(-1i*omegan(k)*x/cn(k)) * ...
            (1.0 - Gamma(k)*exp(1i*omegan(k)*2.0*(x-L)/cn(k)))/Z0(k);
    end
    
end