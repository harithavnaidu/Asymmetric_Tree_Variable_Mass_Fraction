%==========================================================================
% Newton-Raphson method to get mass at each segment
function [Mtotal, R, phi_e, SMCtoCOL,c,tol] = ...
         mass_optimiz_tree_2(p,Q,k,r0)
% Calls functions:
%   - mechanical_properties_PA_4
%   - mass_fracs_2
%   - dWcdx, dWmdx, dWmdx_act, dWedx
% -----------------
% Code by Hamid, 2018 
% Last modified by Vasilina on April 14, 2020
%==========================================================================
%
global rho_w beta gamma R0 alpha_t Ghe1 Ghe2 Ghc Ghm  fiber_ang   

%% Initialization
lhat = r0/R0;

LamZ = 1;

%% Newton-Raphson optimization
err = 1000;
Y = lhat; %

c=0; 
tol=1e-6;
while err>=tol   
    kc = mechanical_properties_PA_4(p); %also updates Ghe1 Ghe2 Ghc
    
    lhat = Y;
    [phi_e, phi_c, phi_m] = mass_fracs_2(2*lhat*R0);
    
    SMCtoCOL = phi_m/phi_c;
    
    beta_t = 1/LamZ*(1/(1+SMCtoCOL)*(0.2+ 0.3*sin(fiber_ang)^2 + ...
            0.3*sin(-fiber_ang)^2)*Ghc*dWcdx(Ghc,kc) +...
            SMCtoCOL/(1+SMCtoCOL)*(Ghm*dWmdx(Ghm,kc)+dWmdx_act(1)));
    
    phi_t = 1 - phi_e;
    
    alpha_act = 0.00872;
    
    VSM_act = dWmdx_act(1)*0.3*rho_w;
    
    alpha_total = alpha_t + alpha_act*VSM_act;
    alpha_act = 0;
        
    beta_e = 1/LamZ*dWedx(1, Ghe1, Ghe2,kc);
    %    active-tone cost +2*pi*alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(phi_e*beta_e+phi_t*beta_t)*(3*p*lhat^2*R0^3)+...
    DF = 2*pi/(0.3*rho_w)*alpha_total*SMCtoCOL/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*lhat*R0^2)+ ...
        2*pi/(0.3*rho_w)*alpha_t*1/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*lhat*R0^2)+...
        alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(0.3*rho_w) + ...
        2*beta*R0^2*lhat - 4*gamma*Q^2*R0^(-4)*lhat^(-5);
    
    D2F = 2*pi/(0.3*rho_w)*alpha_total*SMCtoCOL/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*R0^2)+...
        2*pi/(0.3*rho_w)*alpha_t*1/(SMCtoCOL+1)*phi_t/(phi_e*beta_e+phi_t*beta_t)*(2*p*R0^2)+...
        2*pi*alpha_act*SMCtoCOL/(SMCtoCOL+1)*phi_t*dWmdx_act(Ghm)/(phi_e*beta_e+phi_t*beta_t)*(6*p*lhat*R0^3)+...
        2*beta*R0^2 + 20*gamma*Q^2*R0^(-4)*lhat^(-6);
    
    sol = DF/D2F;
    err = sqrt( sum(sol.^2)/sum(Y.^2) );
    Y = Y - sol;
    c=c+1;
end

lhat = Y;

[phi_e, phi_c, phi_m] = mass_fracs_2(2*lhat*R0);

SMCtoCOL = phi_m/phi_c;
beta_t = 1/LamZ*(1/(1+SMCtoCOL)*(0.2+ 0.3*sin(fiber_ang)^2 ...
       + 0.3*sin(-fiber_ang)^2)*Ghc*dWcdx(Ghc,kc)...
       + SMCtoCOL/(1+SMCtoCOL)*(Ghm*dWmdx(Ghm,kc)+dWmdx_act(1)));
phi_t = 1 - phi_e;

Mtotal = (p*R0*lhat)/(phi_e*beta_e+phi_t*beta_t);
R = R0*lhat;

%VF:
if (R < 0)
    disp('ERROR: Negative radius in optimization. Exit program.');
    exit;
end

% disp(['     Mass-Radius Optimization at vessel ',num2str(k),...
%      ': Newton-Raphson iterations # ',num2str(c),', tollerance ',...
%       num2str(tol),]);
