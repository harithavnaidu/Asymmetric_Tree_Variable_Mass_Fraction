%==========================================================================
% Input impedance of a segment, in frequency domain
function [Zinp] = InputImpedanceOneVessel(NumModes,omegan,L,Z0,c,Gamma)
% only for perturbations
% Zinp(1)=0 is steady part
% ---------
% by Vasilina, 2018
%==========================================================================
    Zinp = zeros(1,NumModes);
    for k=2:NumModes
        Zinp(k) = Z0(k) * (1.0 + Gamma(k)*exp(-1i*omegan(k)*2.0*L/c(k)))...
            /(1.0 - Gamma(k)*exp(-1i*omegan(k)*2.0*L/c(k)));  
    end
   
end