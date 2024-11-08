function [TTs, varargout] = data2TT(Dirname)
% DATA2TT
% Read all the data from Dirname and return it as a cell array of 
% timetables. We furthermore assume that the filename is specified as
% `metadata1_metadata2_metadata3` and also return a cell array.
% Example:
% [TTs, modelname, ref_load_name, turbulence_name] = data2TT(Dirname);

% First get the names of the csv files and outb files.
fnames_csv = dir(fullfile(Dirname, '**/*.csv'));

% For each, get the timetable and naming information.
N = length(fnames_csv);
M = length(split(fnames_csv(1).name, '_'));
TTs = cell(N, 1);
metadata = cell(M, N);
for i = 1:N
    TTs{i} = getTT(fnames_csv(i));
    [~, fname, ~] = fileparts(fnames_csv(i).name);
    metadata(:, i) = split(fname, '_');
end

varargout = cell(M, 1);
for i = 1:M
    varargout(i) = {metadata(i, :)'};
end

end

function TT = getTT(fname_csv)
% getTT
% Read the csv and outb data and combine them into a single timetable.

fname_csv = fullfile(fname_csv.folder, fname_csv.name);
fname_outb = replace(fname_csv, '.csv', '.SFunc.outb');

TT1 = readtimetable(fname_csv);
% We've saved the time data as hh:mm:ss.SSS, but it's nicer to just see the
% seconds on the x-axis, so convert the format to seconds.
TT1.Time.Format = 's';

[Outdata, OutList] = ReadFASTbinary(fname_outb);
TT2 = out2TT(Outdata, OutList);

TT = [TT1, TT2];
end

% function model = getmodel(fname)
% fname_split = split(fname, '_');
% model = fname_split{2};
% end