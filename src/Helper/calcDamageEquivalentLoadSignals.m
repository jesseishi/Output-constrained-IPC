function T = calcDamageEquivalentLoadSignals(TT, signals, m)
% CALCDAMAGEEQUIVALENTLOAD
% TODO: We've got a bunch of similar functions to this. Shouldn't we wrap
% these with one outer function called calcFuncSignals(TT, signals, func)
% so that takes a function (e.g. calcDEL) so that we have less repetitive
% code.
arguments
    TT timetable
    signals string  % Typically a string array.
    m double = 10  % Wöhler slope (typically 4 for steel and 10 for composites).
end

Ttotal = seconds(TT.Time(end) - TT.Time(1));

T = table;
for i = 1:length(signals)
    signal = signals(i);

    c = rainflow(TT(:, signal));
    count = c(:, 1);
    stress_range = c(:, 2);
    DEL = calcDEL(count, stress_range, Ttotal, m);

    newSignalName = sprintf('%s_DEL', signal);
    T.(newSignalName) = DEL;
end
end

function DEL = calcDEL(n, R, neq, m)
DEL = (sum(n .* R.^m) / neq)^(1/m);
end
