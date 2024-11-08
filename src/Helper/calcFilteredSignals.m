function TT = calcFilteredSignals(TT, signals, opts)
%CALCFILTEREDSIGNALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    TT timetable
    signals string  % typically a string array
    opts.signalNameHasUnit logical = false
    opts.unitDelimiter string = "_"
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
    if opts.signalNameHasUnit
        temp = split(signal, opts.unitDelimiter);
        basename = join(temp(1:end-1), opts.unitDelimiter);
        unit = temp(end);
        newSignalName = sprintf('%s_%s_%s', basename, 'filt', unit);
    else
        newSignalName = sprintf('%s_filt', signal);
    end
    
    TT.(newSignalName) = filtfilt(b, a, TT.(signal));
%     TT.(newSignalName) = filtfilt(d, TT.(signal));
end

