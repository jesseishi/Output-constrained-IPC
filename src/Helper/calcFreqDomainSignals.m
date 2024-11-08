function T_fft = calcFreqDomainSignals(TT,signals,opts)
% CALCFREQDOMAINSIGNALS 
% Calclate the frequency domain of each signals in cell array 'signals' of
% the timetable TT and return a new timetable with those new signals.
arguments
    TT timetable
    signals string  % typically a string array
    opts.spectrumtype {mustBeMember(opts.spectrumtype, ["psd", "power"])} = 'psd'
end

% Get the average timestep of the data to use later in the fft.
DT = mean(diff(seconds(TT.Time)));

% Now build a table with the frequency domain signals.
T_fft = table;
for i = 1:length(signals)
    signal = signals(i);
    newSignalName = sprintf('%s_%s', signal, opts.spectrumtype);

    [T_fft.(newSignalName), f] = calcFreqDomainSignal(TT, signal, DT, opts.spectrumtype);
end

% Also add the frequency to the table. Also add the normalized frequency,
% so that you can easily see the 1P/2P, etc...
f_1P_rpm = mean(TT.GenSpeed);
f_1P_hz = rpm2hz(f_1P_rpm);
T_fft.f_hz = f;
T_fft.f_1P = f ./ f_1P_hz;
end

function [Pxx, f] = calcFreqDomainSignal(TT, signal, DT, spectrumtype)

    x = TT.(signal);
    x = detrend(x, 'constant');

    % TODO: windowing.
    % When specifying the spectrumtype to be 'psd' the resulting unit it
    % signalunit^2/Hz and when the spectrumtype is 'power' the resulting
    % unit is signalunit^2.
    [Pxx, f] = periodogram(x, [], [], 1/DT, spectrumtype);
end
