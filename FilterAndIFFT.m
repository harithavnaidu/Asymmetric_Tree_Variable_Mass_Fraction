%==========================================================================
% Get function transformation from frequency to time domain: F -> f
function [fTotal,fOscil]=FilterAndIFFT(NumModes,N,Fn)
% Use filtering of messy modes(need to check for each flow):
%    - clean from zero Fn, excluding steady.
%    - truncate 10 mods.
% N > 2*NumModes-2
% Total = steady state + oscillatory contribution
% ---------
% by Vasilina, 2018
%==========================================================================
    F = zeros(N,1);
    F(1:NumModes) = Fn;
    
    % symmetric conjugate
    % F(k) = F(N-k+2), k=2,...N, F(2)=F(N), F(10)=F(N-8)
    F(N-NumModes+2:N)=flip(conj(Fn(2:NumModes)));
 
    fTotal = N*ifft(F,'symmetric'); 
    fOscil = fTotal - F(1);
end