%==========================================================================
% This function uses SoL resutls to find the Young's modulus of an 
% inflated thin-walled cylinder, given the masses.
% l1 and l2 are 1, unless another point is needed for linearizaiton.
%
function [Ett, Ezz, nu_tz, nu_zt, C] = ...
         YoungMod_2(Me, Mt, lt, lz,p,SMCtoCOL)
% Calls function:
%   - mechanical_properties_PA_4
% -----------------
% Code by Hamid, 2018 
% Format updates by Vasilina, April 16, 2020
% =========================================================================
%
global Ghe1 Ghe2 Ghc Ghm lmin lmax fiber_ang S rho_w

fiber_angle = [0, pi/2, -fiber_ang, fiber_ang];
l2 = lz; l1 = lt; % The reference configuration is the homeostatic condition (F=I)
% lkc2 = Ghc^2;%Ghc^2*(l2^2*cos(fiber_angle).^2 + l1^2*sin(fiber_angle).^2);
lkm2 = Ghm^2;

% Material Properties
kc = mechanical_properties_PA_4(p);

%%%%%
Kact = S/rho_w;
%%%%%

ddwcddx = [0,0,0,0];exp_Q = [0,0,0,0];lkc2 = [0,0,0,0];
for i = 1:4
    lkc2(i) = Ghc^2*(l1^2*sin(fiber_angle(i))^2 + l2^2*cos(fiber_angle(i))^2);  
    exp_Q(i) = exp(kc(3)*(lkc2(i)-1.0).^2);
    ddwcddx(i) = 1/2*kc(2)*(2*kc(3)*(lkc2(i)-1).^2+1).*exp_Q(i);
end

exp_M = exp(kc(5)*(lkm2-1.0).^2);
ddwmddx = 1/2*kc(4)*(2*kc(5)*(lkm2 - 1).^2+1).*exp_M;
ddTddx = (Kact*(-lt + lmax))/(lt*(lmax - lmin)^2);


Ms = SMCtoCOL/(SMCtoCOL+1)*Mt;
Mc = 1/(SMCtoCOL+1)*Mt;
Mk = [0.2, .2, 0.3, 0.3]*Mc; %% guessed!

ACT = ddTddx;
Atttt = sum(Ghc^4*Mk.*ddwcddx.*sin(fiber_angle).^4) + Ghm^4*Ms*ddwmddx + Ms*ACT;
Attzz = sum(Ghc^4*Mk.*ddwcddx.*(sin(fiber_angle).^2).*cos(fiber_angle).^2);
Azzzz = sum(Ghc^4*Mk.*ddwcddx.*cos(fiber_angle).^4);

Tt = Me*kc(1)*Ghe1^2*lt^2 + kc(2)*lt^2*sum(Mk.*(lkc2.*(lkc2 - 1).*exp(kc(3)* (lkc2 - 1).^2).*sin(fiber_angle).^2)) + ...
    Ms*lt^2*kc(4)*lkm2*(lkm2 - 1)*exp(kc(5)*(lkm2 - 1)^2) +Ms*lz*(Kact*(1 - (lmax - lt)^2/(lmax - lmin)^2));% + 2*kjterm;%)/((Mt+Me)/rho_w

Tz = Me*kc(1)*Ghe2^2*lz^2 + lz^2*sum(Mk.* (kc(2)*lkc2.*(lkc2 - 1).*exp(kc(3)* (lkc2 - 1).^2).*cos(fiber_angle).^2));% + 2*kjterm;%)/((Mt+Me)/rho_w

Ctttt = Tt + 2*Atttt*lt^4;
Czzzz = Tz + 2*Azzzz*lz^4;
Cttzz = 2*Attzz*lt^2*lz^2;

C = [Ctttt,Cttzz;Cttzz,Czzzz]./((Mt+Me)/(0.3*rho_w));

Sttzz = -Cttzz/(Ctttt*Czzzz - Cttzz^2)*((Mt+Me)/(0.3*rho_w));

Ett = (Ctttt*Czzzz - Cttzz^2)/Czzzz/((Mt+Me)/(0.3*rho_w));
Ezz = (Ctttt*Czzzz - Cttzz^2)/Ctttt/((Mt+Me)/(0.3*rho_w));

nu_tz = -Ett*Sttzz;
nu_zt = -Ezz*Sttzz;
end
