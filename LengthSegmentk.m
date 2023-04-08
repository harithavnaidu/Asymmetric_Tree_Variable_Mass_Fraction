%==========================================================================
% Get branch length as a function of radius
function [L]=LengthSegmentk(R)
% ------
% by Vasilina, 2018
%==========================================================================
    % from [Olufsen 2012] for 4-12 order in Huang paper
    % empirical relations given in mm = 0.001m
    LOluf = 12.4*(R*1000).^1.1/1000;
    
    % scale by factor half to get more physiological length (shorter)
    % as seen in larger vessels
    L = LOluf/2; 
end