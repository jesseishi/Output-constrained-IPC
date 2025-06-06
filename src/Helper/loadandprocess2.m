function [TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(nameFilter,opts)
%LOADANDPROCESS2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    nameFilter string  % Can also be a string array to combine multiple searches
    opts.TStart duration = seconds(100)
    opts.TEnd duration = seconds(Inf)
    opts.verbose logical = true
    opts.saveTTs logical = true
    opts.saveTsfft logical = true
    opts.saveT logical = true
    opts.assertEqualLength logical = true
    opts.signalNotFound = 'none'
end

fnames = getFnames(nameFilter);
S = timerange(opts.TStart, opts.TEnd);

% Loop through the files and extract the necessary data.
N = length(fnames);
M = length(split(fnames(1).name, '_'));
metadata = cell(M, N);
TTs = cell(N, 1);
TMaxs = cell(N, 1);
Ts_fft = cell(N, 1);
Ts = cell(N, 1);
% TODO: use parfor to make this faster. + use Matlab's profiler to check
% where this code is slow.
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool('Processes')
end

parfor i = 1:N
    fname = fnames(i);

    if opts.verbose
        fprintf("%02d: Loading %s\n", i, fname.name);
    end
%     continue

    % Let's start by extracting metadata from the filename.
    [~, name, ~] = fileparts(fname.name);
    metadata(:, i) = split(name, '_');

    % Now we load the timetable, store the maximum time in TMaxs, and do
    % some postprocessing.
    TT = getTT(fname);
    TT = calcConvertedSignals(TT, 'deleteOldName', true, 'signalNotFound', opts.signalNotFound);
    TT = calcMagnitudeAndPhaseSignalsPairs(TT, ["M_tilt_MNm", "M_yaw_MNm", "theta_tilt_deg", "theta_yaw_deg", "M_tilt_0_MNm", "M_yaw_0_MNm"], ["M_ty", "theta_ty", "M_ty_0"], 'signalNotFound', opts.signalNotFound, 'signalNameHasUnit', true);
    TT = calcFilteredSignals(TT, ["M_ty_mag_MNm", "theta_ty_mag_deg", "M_tilt_MNm", "M_yaw_MNm", "theta_tilt_deg", "theta_yaw_deg"], 'signalNameHasUnit', 'true');
    TT = calc1PFilteredSignals(TT, ["FlappingMoment1_MNm", "BldPitch1_deg"], 'signalNameHasUnit', true, 'TStart', seconds(100));

    % Let's now calculate some results in the frequency domain.
%     T_fft = calcFreqDomainSignals(TT(S,:), ["FlappingMoment1_MNm", "BldPitch1_deg"], 'spectrumtype', 'psd');

    % And finally we calculate some metrics for this entire run.
    Tnew = [calcDutyCycleSignals(TT(S,:), ["BldPitch1_deg", "BldPitch2_deg", "BldPitch3_deg"], 2), ...
            calcDamageEquivalentLoadSignals(TT(S,:), ["FlappingMoment1_MNm", "FlappingMoment2_MNm", "FlappingMoment3_MNm"]), ...
            ... calc1PPowerSignals(T_fft, ["FlappingMoment1_MNm_psd", "BldPitch1_deg_psd"]), ...
            calcMeanSignals(TT(S,:), ["theta_ty_mag_deg", "M_ty_mag_MNm"])];

    % Store the data if we want to. We might not want to save all
    % timetables and only look at final metrics when we are looking at a
    % lot of data.
    TMaxs{i} = TT.Time(end);
%     if opts.assertEqualLength
%         assert(TMaxs{1} == TMaxs{i}, "All simulations must have the same length")
%     end
    if ~isinf(opts.TEnd) && (TMaxs{i} < opts.TEnd)
        warning('Option TEnd set to %f but only have %f of data', seconds(opts.TEnd), seconds(TMaxs{i}))
    end
    if opts.saveTTs
        TTs{i} = TT;
    end
%     if opts.saveTsfft
%         Ts_fft{i} = T_fft;
%     end
    if opts.saveT
        Ts{i} = Tnew;
    end
end

T = vertcat(Ts{:});  % Join all the individual tables together.

end


% Function to get all filenames from multiple name filters.
function fnames = getFnames(nameFilter)
arguments
    nameFilter string
end

% First get the names of the csv files.
fnames = dir(nameFilter(1));
for iFile = 2:length(nameFilter)  % Looping through 2:N will skip if N == 1.
    temp = dir(nameFilter(iFile));
    fnames = [fnames; temp];
end
% fnames_csv = dir(nameFilter);
if isempty(fnames)
    error("Couldn't find any csv's in %s", nameFilter)
end
end


% Function to read a csv file to a timetable object.
function TT = getTT(fname_csv)
fname_csv = fullfile(fname_csv.folder, fname_csv.name);
TT = readtimetable(fname_csv);
% We've saved the time data as hh:mm:ss.SSS, but it's nicer to just see the
% seconds on the x-axis, so convert the format to seconds.
TT.Time.Format = 's';
end

