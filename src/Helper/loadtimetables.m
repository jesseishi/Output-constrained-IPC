function [TTs, N, metadata] = loadtimetables(searchString, x)
% loadtimetables
% Read all the data from Dirname and return it as a cell array of 
% timetables. We furthermore assume that the filename is specified as
% `metadata1_metadata2_metadata3` and also return a cell array.
% Example cIPCBox_ref200_turb10.csv
% [TTs, modelname, ref_load_name, turbulence_name] = data2TT(Dirname);
arguments
    searchString string  % Can be a single string or a string array (in which case multiple results will be combined).
    x.verbose logical = true
end

% First get the names of the csv files.
fnames_csv = dir(searchString(1));
for i = 2:length(searchString)  % Looping through 2:N will skip if N == 1.
    temp = dir(searchString(i));
    fnames_csv = [fnames_csv; temp];
end
% fnames_csv = dir(nameFilter);
if isempty(fnames_csv)
    error("Couldn't find any csv's in %s", searchString)
end

% For each, get the timetable and naming information.
N = length(fnames_csv);
M = length(split(fnames_csv(1).name, '_'));
TTs = cell(N, 1);
metadata = cell(M, N);
for i = 1:N
    if x.verbose
        fprintf("Loading %s\n", fnames_csv(i).name);
    end
    TTs{i} = getTT(fnames_csv(i));
    [~, fname, ~] = fileparts(fnames_csv(i).name);
    metadata(:, i) = split(fname, '_');
end

% varargout = cell(M, 1);
% for i = 1:M
%     varargout(i) = {metadata(i, :)'};
% end

end

function TT = getTT(fname_csv)
fname_csv = fullfile(fname_csv.folder, fname_csv.name);
TT = readtimetable(fname_csv);
% We've saved the time data as hh:mm:ss.SSS, but it's nicer to just see the
% seconds on the x-axis, so convert the format to seconds.
TT.Time.Format = 's';
end
