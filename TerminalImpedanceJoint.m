%==========================================================================
% Get terminal impedance of parent vessel computed at joint
function [ZT] =  TerminalImpedanceJoint(NumModes,Zinp1,Zinp2)
% ---------
% by Vasilina, 2018
%==========================================================================  
    ZT = zeros(1,NumModes);
    for k = 2:NumModes
        if ((Zinp1(k) + Zinp2(k))==0.0 )
            disp('exit: zero = Zinp1+Zinp2');
            break
        else
            ZT(k)= Zinp1(k)*Zinp2(k)/(Zinp1(k) + Zinp2(k));  
        end
    end    
end