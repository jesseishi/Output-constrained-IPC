function T = calcMeanSignals(TT, signals)
% CALCMEANSIGNALS
arguments
    TT timetable
    signals string
end

T = table;
for i = 1:length(signals)
    signal = signals(i);
    newSignalName = sprintf('%s_mean', signal);
    
    T.(newSignalName) = mean(TT{:, signal});
end

end

