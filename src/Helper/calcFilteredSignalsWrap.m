function TT = calcFilteredSignalsWrap(TT, signals)
%CALCFILTEREDSIGNALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    TT timetable
    signals string  % typically a string array
end

% Numerator (b) and denominator (a) coefficients of transfer function.
w = 1;  % cut-off frequency.
b = [w];
a = [1, w];
sys = tf(b, a);
sysd = c2d(sys, 0.005);
b = sysd.Numerator{:};
a = sysd.Denominator{:};
% TODO: Investigate using designfilt or how to directly design dicrete
% low-pass filters.

for i = 1:length(signals)
    signal = signals(i);
    newSignalName = sprintf('%s_filt', signal);
    
    TT.(newSignalName) = filtfilt(b, a, TT.(signal));
    %     TT.(newSignalName) = filtfilt(d, TT.(signal));
end
end

function y = wrap(x)
y = mod(x + pi, 2*pi) - pi;
end

