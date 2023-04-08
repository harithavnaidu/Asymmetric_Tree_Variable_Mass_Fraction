%==========================================================================
% Function for variable mass fractions of elastin, collagen and SMC
function [e, c, m] = mass_fracs_2(D)
% -----------------
% Code by Hamid, June 13, 2018 
%==========================================================================
% Elastin

% no transitin, all elastic arteries
%  x1 = 100;
%  x2 = 90;
% 
x1 = 4200; %3200;
x2 = 1000; %2000;

y1 = 11;
y2 = 2;

% coefficients = polyfit([x1, x2], [y1, y2], 1);
% ae = coefficients (1);
ae = (y1-y2)/(x1-x2);
be = y1+y2;

% SMCs

y3 = 12;
y4 = 39;

% coefficients = polyfit([x1, x2], [y3, y4], 1);
% am = coefficients (1);
am = (y3-y4)/(x1-x2);
bm = y3+y4;

D = D*10^6;
e = elastin(D,ae,y1,y2,x1,x2)/100;
m = smc(D,am,y3,y4,x1,x2)/100;
c = 1 - e - m;


function e = elastin(x,ae,y1,y2,x1,x2)
    e = y2+(y1-y2)*(atan(ae*(x-(x1+x2)/2))+pi/2)/pi;
end

function m = smc(x,am,y3,y4,x1,x2)
    m = y4+(y3-y4)*(atan(-am*(x-(x1+x2)/2))+pi/2)/pi;
end

end