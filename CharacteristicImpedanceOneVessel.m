%==========================================================================
% Womersley solution for characteristic impedance
function [Zn] = CharacteristicImpedanceOneVessel(NumModes,R,rho,cn,Mn,gn)
% get complex values in frequency domain
% ---------
% by Vasilina, 2018
%==========================================================================
    Zn = zeros(1,NumModes); 
    for k = 2:NumModes
        Zn(k)=rho*cn(k)/(pi*R^2*(1.0-Mn(k)*gn(k)));
    end
end