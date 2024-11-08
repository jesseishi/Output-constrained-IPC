function T = calcDutyCycleSignals(TT, signals, signals_dot_max)
% CALCDUTYCYCLE

arguments
    TT timetable
    signals string  % Typically a string array.
    signals_dot_max double = 1;
end

DT = mean(diff(seconds(TT.Time)));
if size(signals_dot_max) == 1
    signals_dot_max = signals_dot_max * ones(size(signals));
elseif size(signals_dot_max) ~= size(signals)
    error('Size of `signals_dot_max` must be 1 or match the size of `signals`.')
end

T = table;
for i = 1:length(signals)
    signal = signals(i);
    signal_dot_max = signals_dot_max(i);

    DC = calcDutyCycleSignal(TT{:, signal}, DT, signal_dot_max);
    newSignalName = sprintf('%s_DC', signal);
    T.(newSignalName) = DC;
end
end


function ADC = calcDutyCycleSignal(x, DT, x_dot_max)
% CALCACTUATORDUTYCYCLE
arguments
    x double
    DT
    x_dot_max double
end
Ttot = DT * length(x);
x_dot = abs(diff(x)) ./ DT;
ADC = 1/Ttot * trapz(x_dot) * DT / x_dot_max;
end
