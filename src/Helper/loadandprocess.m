function [TTs, Ts_fft, T, S, N, metadata] = loadandprocess(nameFilter, opts)
arguments
    nameFilter string  % Can be a single string or a string array (in which case multiple results will be combined).
    opts.TStart duration = seconds(100)
end

% Load all timetables from nameFilter.
[TTs, N, metadata] = loadtimetables(nameFilter);

% Do some postprocessing.
TMax = TTs{1}.Time(end);
for i = 2:N
    assert(TMax == TTs{i}.Time(end), 'Not all simulations are the same length')
end
S = timerange(opts.TStart, TMax);

% Calculate some derived signals and add them to the timetable.
for i = 1:N
    TTs{i} = calcConvertedSignals(TTs{i}, 'deleteOldName', true, 'signalNotFound', 'warning');
    TTs{i} = calcMagnitudeAndPhaseSignalsPairs(TTs{i}, ["M_tilt_kNm", "M_yaw_kNm", "theta_tilt_deg", "theta_yaw_deg", "M_tilt_0_kNm", "M_yaw_0_kNm"], ["M_ty", "theta_ty", "M_ty_0"], 'signalNotFound', 'warning', 'signalNameHasUnit', true);
    TTs{i} = calcFilteredSignals(TTs{i}, ["M_ty_mag_kNm", "theta_ty_mag_deg", "M_tilt_kNm", "M_yaw_kNm", "theta_tilt_deg", "theta_yaw_deg"], 'signalNameHasUnit', 'true');
end

% Calculate some results in the frequency domain.
Ts_fft = cell(N);
for i = 1:N
    Ts_fft{i} = calcFreqDomainSignals(TTs{i}, ["FlappingMoment1_MNm", "BldPitch1_deg"], 't_steadystate', opts.TStart, 'spectrumtype', 'psd');
end

% Calculate some metrics for each run.
T = table;
for i = 1:N

    % TODO: Combine into a single function? These now have a ton of overlap
    % between them.
    T = [T;
         calcDutyCycleSignals(TTs{i}, ["BldPitch1_deg", "BldPitch2_deg", "BldPitch3_deg"], opts.TStart, 1.0), ...
         calcDamageEquivalentLoadSignals(TTs{i}, ["FlappingMoment1_MNm", "FlappingMoment2_MNm", "FlappingMoment3_MNm"], opts.TStart), ...
         calc1PPowerSignals(Ts_fft{i}, ["FlappingMoment1_MNm_psd", "BldPitch1_deg_psd"]), ...
         calcMeanSignals(TTs{i}, ["theta_ty_mag_deg", "M_ty_mag_kNm"], opts.TStart)];
end

%% Check that all results have the same length.
N1 = length(TTs);
N2 = length(Ts_fft);
N3 = height(T);
assert(N1 == N2)
assert(N2 == N3)
N = N1;


end


