%**************************************************************************
% Function to find unstressed radius
% - ZeroP
% Internal functions:
% - ddWeddLn
% - ddWkddLn
% - dWedLn
% - dWkdLn
% - dWmdLn
% - ddWmddLn
%--------------------------------------------------------------------------
% Code by Hamid, 2018
% Updated format by Vasilina, April 16, 2020
%**************************************************************************
%% =========================================================================
function r0 = ZeroP(Me,Mt,R,Pmid,P)
% Calls external functions:
%   - mass_fracs_2
%   - mechanical_properties_PA_4
%   - dWmdx_act
%==========================================================================
global Ghe1 Ghe2 Ghc Ghm fiber_ang 

[phi_e, phi_c, phi_m] = mass_fracs_2(2*R);
SMC2COL = phi_m/phi_c;


Mm = SMC2COL/(SMC2COL+1)*Mt;
Mc = 1/(SMC2COL+1)*Mt;
Mk = [0.2, 0.2, 0.3, 0.3]*Mc;

Lz = 1.0;
% P = 0;
% Material Properties
kc = mechanical_properties_PA_4(Pmid); %also updates Ghe1 Ghe2 Ghc
%
r0 = R;
err = 1000;

% NR iterations
while err > 10^-6
    Lt = r0/R;
    % h = Mtotal/(rho_w*Lt*Lz);
    
    % elastin: isotropic matrix
    Lt_e = Ghe1*Lt;
    Lz_e = Ghe2*Lz;
    
    dwdLt = Me * Ghe1 * dWedLn(Lt_e,Lz_e, 1,kc);
    d2wddLt = Me * Ghe1^2 * ddWeddLn(Lt_e,Lz_e, 1,kc);
    
    % collagen: four fiber families
    alpha = [0, pi/2, atan(Lt/Lz*tan(fiber_ang)), -atan(Lt/Lz*tan(fiber_ang))];
    dwdLtc = 0;
    d2wddLtc = 0;
    for k = 1:4
        
        Lk = sqrt(Ghc^2*(Lt^2*sin(alpha(k))^2+Lz^2*cos(alpha(k))^2));
        
        dwdLtc = Mk(k) * ( Ghc^2 * Lt * sin(alpha(k))^2/Lk ) * dWkdLn(Lk,kc) + dwdLtc;
        
        d2wddLtc = Mk(k) * ( (2* Ghc^2 * Lt * sin(alpha(k))^2/Lk )^2 * ddWkddLn(Lk,kc) + ...
            Ghc^4 * ( 2*Lz^2*sin(alpha(k))^2*cos(alpha(k))^2/Lk^3 ) *  dWkdLn(Lk,kc) ) + d2wddLtc;
        
    end
    
    dwdLt = dwdLt + dwdLtc;
    d2wddLt = d2wddLt + d2wddLtc;
    
    % SMC: oriented circumferentially, passive and active response
    % Passive response
    Lm = Ghm*Lt;
    
    dwdLtm = (Mm) * ( Ghm) * dWmdLn(Lm,kc);
    
    d2wddLtm = (Mm) * (2 * Ghm^2 * ddWmddLn(Lm,kc));
    
    dwdLt = dwdLt + dwdLtm;
    d2wddLt = d2wddLt + d2wddLtm;
    
    % Active muscle tone (Constant)
    T_act = Mm*dWmdx_act(1);
    
    % Function and derivaite
    F = 1/Lz * dwdLt + T_act - P * r0;
    dFdr = 1/R * 1/Lz * d2wddLt - P;
    
    r_new = r0 - F/dFdr;
    err = abs((r_new-r0)/r0);
    r0 = r_new;
    
end

%% ======================================================================== 
% Derivatives of strain energy function
    function  y=ddWeddLn(Ln_t, Ln_z, f, kc)
%=========================================================================      
        if f==1
            y = kc(1) * (1.0+3.0/(Ln_t^4*Ln_z^2));
        elseif f==2
            y = kc(1) *(1.0+3.0/(Ln_z^4*Ln_t^2));
        elseif  f==3
            y = kc(1) * 2.0 /(Ln_t*Ln_z)^3;
        else
            exit('Wrong parameter for ddWddLn');
        end
        
    end

%% ========================================================================
    function  y=ddWkddLn(Ln,kc)
% ========================================================================     
        exp_Q = exp(kc(3)*(Ln^2-1.0)^2);
        y = kc(2)* (3*Ln^2-1.0+4*kc(3)*(Ln*(Ln^2-1.0))^2)*exp_Q;
        
    end

%% ========================================================================
    function  y=ddWmddLn(Ln,kc)
% =========================================================================    
        exp_Q = exp (kc(5)*(Ln^2-1.0)^2);
        y=kc(4)*(3*Ln^2-1.0+4*kc(5)*(Ln*(Ln^2-1.0))^2)*exp_Q;
    end

%% ========================================================================
     function  y=dWedLn(Ln_t, Ln_z, axis,kc)
% =========================================================================
        if axis == 1  %circumferential
            y=kc(1)*(Ln_t - 1.0/(Ln_t^3*Ln_z^2));
        elseif axis ==2  %axial
            y= kc(1)*(Ln_z - 1.0/(Ln_z^3*Ln_t^2));
        else
            exit('wrong parameter in dWedLn');
        end
    end

%% ========================================================================
    function  y=dWkdLn(Ln,kc)
% =========================================================================    
        y=kc(2)*Ln*(Ln^2-1.0)*exp(kc(3)*(Ln^2-1.0)^2);
    end
%% ========================================================================
    function  y=dWmdLn(Ln,kc)
% =========================================================================  
        y= kc(4)*Ln*(Ln^2-1.0)*exp(kc(5)*(Ln^2-1.0)^2);
    end
end