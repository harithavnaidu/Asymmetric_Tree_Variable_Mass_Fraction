%==========================================================================
% Compute mechanical properties for each segment based on mid-pressure
function KC = mechanical_properties_PA_4(p)
% Get: KC = [c1, c2, c3, c4, c5]
% NOTE: It also updates prestreches Ghe1 Ghe2 Ghc 
%       which are calibrated w.r.t. two experimental points in h/D curve.
%       P_h is a parameter, set to specify homeostatic pressure at loaded 
%       state. Coeff. and exponent for 'corr' are calibrated.
% -----------------
% Code by Hamid, 2018 
%   - Updated by Vasilina, January 13, 2021
%==========================================================================
%
 global Ghe1 Ghe2 Ghc Ghm rho_w P_h ...
        %pruneflag  PTflag
 
% Passive wall

% correction factor for changes of mid-pressure w.r.t homeostatic
corr = p/P_h;

%elastin prestretches
Ghe1 = 1.27*corr^0.15;   % nondim.
% if (Ghe1 < 1) %compression
%     Ghe1 = 1.2; %VF
% end
Ghe2 = Ghe1;            % nondim.

% from experiments (Wang,Hamid... 2020)
c1 = 28.8341;           % for elastin
c2 = 178.5934;          % for collagen
c3 = 1.05;              % for collagen
c4 = 24.5093;           % for passive smooth muscle
c5 = 0.75;              % for passive smooth muscle
% c1 =  Sigma_e / ((Ghe1^2-(Ghe1^2*Ghe2^2)^(-1))*(0.3*rho_w));

% fiber_ang = 45*pi/180; % angle of helical collagen fibers

% collagen  prestretches
Ghc = 1.154*corr^.15;    % nondim.
% if (Ghc < 1) %compression
%     Ghc = 1.05; %VF
% end

% % for symmetric case verification
% if (pruneflag==0)&&(PTflag==1)
%     Ghc = 1.154*corr^.22;
% end

str_ch = 0.0;
ang = [0.0 90.0 45.0 90+45.0]*(pi/180);
c_frac0 = [.2 .2 .3 .3];
for i=1:4
    str_ch = str_ch + c_frac0(i)*Ghc^2*(Ghc^2-1)*...
        exp(c3*(Ghc^2-1.0)^2)*(sin(ang(i)))^2;
end
% c2 = Sigma_c/ (str_ch*(0.3*rho_w));

% smooth muscle cell prestretch
Ghm = 1.21;             % nondim.

VSM_act = dWmdx_act(1)*0.3*rho_w;
% c4 = (Sigma_m-VSM_act)/(Ghm^2*(Ghm^2-1.0)*exp(c5*(Ghm^2-1.0)^2)*(0.3*rho_w));

% strain energy parameters
KC = [c1, c2, c3, c4, c5];

end