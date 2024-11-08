function TT = calc1PFilteredSignals(TT, signals, opts)
%CALC1PFILTEREDSIGNALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    TT timetable
    signals string  % typically a string array
    opts.signalNameHasUnit logical = false
    opts.unitDelimiter string = "_"
    opts.TStart duration = seconds(0);
end

fs = 1/mean(diff(seconds(TT.Time)));
S = timerange(opts.TStart, seconds(Inf));
onePfreq = rpm2hz(mean(TT.RotSpeed(S)));
fpass = [0.75 * onePfreq, 1.25 * onePfreq];

for i = 1:length(signals)
    signal = signals(i);

    if opts.signalNameHasUnit
        temp = split(signal, opts.unitDelimiter);
        basename = join(temp(1:end-1), opts.unitDelimiter);
        unit = temp(end);
        newSignalName = sprintf('%s_%s_%s', basename, '1Pfilt', unit);
    else
        newSignalName = sprintf('%s_1Pfilt', signal);
    end

    % iir uses filtfilt for zero-phase filtering.
    % The beginning of the signal has lots of transients that mess up this
    % filtering so cut it out.
    yFilt = bandpass(TT{S, (signal)}, fpass, fs, 'ImpulseResponse', 'iir', 'Steepness', 0.95, 'StopBandAttenuation', 60);
    yFilt = [zeros(seconds(opts.TStart) * fs, 1); yFilt];

    % That's also not nice. Let's investigate again.
    yFilt = bandpass(TT{:, (signal)}, fpass, fs);
    TT.(newSignalName) = yFilt;
end

end