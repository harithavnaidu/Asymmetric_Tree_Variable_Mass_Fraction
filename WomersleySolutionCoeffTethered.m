%==========================================================================
% Get Womersley solution coefficients for longitudinally constrained
% elastic wall (Womersley 1957)
function [cn,Mn,gn,c_Rn,alphan] = WomersleySolutionCoeffTethered...
                                (NumModes,omegan,mu,rho,R,a,h)
% Get fequency-dependent and complex values of wave speed and M,g coeff.        
% Code by Vasilina, July 2018
%   - Updated April 15, 2020
%==========================================================================
cn = zeros(NumModes,1);
Mn = zeros(NumModes,1);
gn = zeros(NumModes,1);

% first component is dummy
% Inviscid wave speed: c_0 = sqrt(E*h/(2*rho*R));

%frequency, Womersley num, Lambda and g-factor vectors
alphan = zeros(NumModes,1);
Lambdan = zeros(NumModes,1);
for k=2:NumModes
    alphan(k) = sqrt(rho*omegan(k)/mu)*R;
    Lambdan(k,1) = alphan(k,1)*1i^1.5;
    gn(k) = 2.0*besselj(1,Lambdan(k))/(Lambdan(k)*besselj(0,Lambdan(k)));
end

% Wave speeds complex solution with stiffness component (a)
for k=2:NumModes
%   cn(k) = c_0*sqrt((1.0-gn(k))/(1.0-sigma^2)); 
    cn(k) = sqrt(a*h*(1.0-gn(k))/(2*rho*R)); 
end

% Wave speed real and imaginary parts
inverse_cn = zeros(NumModes-1,1);
c_Rn = zeros(NumModes,1);
c_In = zeros(NumModes,1);
for k=2:NumModes
    inverse_cn(k) = cn(k)^(-1);
    c_Rn(k) = 1.0/real(inverse_cn(k));
    c_In(k) = 1.0/imag(inverse_cn(k));
end
% Estimate of the spatial wavelengt: Lwave = c_Rn(1)*T;

% Get M (1 for tethered case)
for k=2:NumModes
   Mn(k) = 1.0;
end
end