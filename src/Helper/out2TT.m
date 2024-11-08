function TT = out2TT(OutData,OutList)
% out2TT 
% Out to timetable collection.
[~, n_vars] = size(OutData);

if n_vars ~= length(OutList)
    error('OutList must have the same length as the amount of variables in OutData (assumed to be the 2nd dimension)');
end

TT = array2timetable(OutData(:, 2:end), 'TimeStep', seconds(0.005), 'VariableNames', OutList(2:end));

end
