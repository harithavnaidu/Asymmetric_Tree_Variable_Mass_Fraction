% =========================================================================
% Get function transformation from the time domain to fequency domain
function [Qn] = FilterAndFFT(q,NumModes,N)
%==========================================================================
    Fn = fft(q(1:N)); 
    Qn = zeros(NumModes,1);

    Qn(1) = Fn(1)/N; 
    for k=2:NumModes
        Qn(k) = Fn(k)/N;   
    end
end