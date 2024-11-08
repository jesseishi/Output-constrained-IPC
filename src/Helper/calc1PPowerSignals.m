function T = calc1PPowerSignals(T_fft,signals)
%CALC1PPOWER Summary of this function goes here
%   Detailed explanation goes here
arguments
    T_fft table
    signals string
end

i_1P = find(T_fft.f_1P > 1, 1, 'first');  % Find the first index that's bigger than one.

T = table;
for i = 1:length(signals)
    signal = signals(i);

    newSignalName = sprintf('%s_at1P', signal);
    T.(newSignalName) = calc1PPowerSignal(T_fft.(signal), i_1P);
end
end

function P = calc1PPowerSignal(Pxx, f_approx)

% Get the peak power around this index.
P = max(Pxx(f_approx-5:f_approx+5));

end


