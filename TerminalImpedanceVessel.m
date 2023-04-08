%==========================================================================
% Terminal impedance from reflection coefficient and char. impedance
function [ZT] = TerminalImpedanceVessel(NumModes,Z0,Gamma)
% ---------
% by Vasilina, 2018
%==========================================================================    
    ZT = zeros(1,NumModes);
    for k=2:NumModes
        if (Gamma == 1.0)
            ZT(k)= NaN;
        dips('Warning: vessel has a closed end, ifinite terminal impedance')
        else
            ZT(k) = Z0(k)*(1.0 + Gamma(k))/(1.0 - Gamma(k));
        end
    end    
end