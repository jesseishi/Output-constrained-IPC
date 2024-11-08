function TT = calcMagnitudeAndPhaseSignalsPairs(TT, signal_pairs, signal_basenames, x)
%CALCMAGNITUDEANDPHASESIGNALS Summary of this function goes here
%   Detailed explanation goes here
arguments
    TT timetable
    signal_pairs string  % typically a string array
    signal_basenames string  % typically a string array
    x.signalNotFound string = "error"
    x.signalNameHasUnit logical = true
    x.unitDelimiter string = "_"
end

assert(length(signal_pairs) == 2*length(signal_basenames), "Must have twice as many signal pairs (%i) as basenames (%i).", length(signal_pairs), length(signal_basenames));


for i = 1:length(signal_basenames)
    basename = signal_basenames(i);
    s1 = signal_pairs(2*i-1);
    s2 = signal_pairs(2*i);

    if x.signalNameHasUnit
        temp = split(s1, x.unitDelimiter);
        unit = temp(end);
        temp = split(s2, x.unitDelimiter);
        unit2 = temp(end);
        assert(unit == unit2, "Signal pairs must have the same unit.");

        basename_mag = basename + "_mag_" + unit;
        basename_phase = basename + "_phase_deg";
    else
        basename_mag = basename + "_mag";
        basename_phase = basename + "_phase_deg";
    end


    % Check if signals exist
    if ~all(ismember({char(s1), char(s2)}, fieldnames(TT)))
        switch x.signalNotFound
            case "none"
                continue
            case "warning"
                warning("Couldn't find signals %s %s in timetable.", s1, s2);
                continue
            otherwise
                error("Couldn't find signals %s %s in timetable.", s1, s2);
        end
    end
    
    TT.(basename_mag) = sqrt(TT.(s1).^2 + TT.(s2).^2);
    TT.(basename_phase) = rad2deg(atan2(TT.(s1), TT.(s2)));

end

end
