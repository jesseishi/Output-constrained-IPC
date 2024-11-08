function plotTTsOverview(TTs, opts)
%PLOTTTSOVERVIEW Plot some opts.signals of multiple timetables to check if
%simulations ran correctly.
% Makes N-by-1 subplots where each row has a signal from opts.signals. For
% each subplots the signal from all timetables will be plotted.
arguments
    TTs cell
    opts.legendLabels string = []
    opts.signals string = ["Wind1VelX", "RtSpeed_rpm", "GenTq_kNm", "GenPwr_MW", "theta_CPC_deg", "FlappingMoment1_MNm", "BldPitch1_deg", "M_tilt_kNm", "FlappingMoment1_1Pfilt_MNm"];
    opts.filename string = []
end

% opts.signals = ["RtSpeed_rpm", "GenTq_kNm", "GenPwr_MW", "theta_CPC_deg", "FlappingMoment1_MNm", "TipDyc1_m", "theta_1_deg", "BldPitch1_deg", "M_tilt_kNm", "M_yaw_kNm", ];
[f, axs] = preplot(length(opts.signals), 1, 'interpreter', 'none');
for i = 1:length(opts.signals)
    for j = 1:length(TTs)
        plot(axs(i), TTs{j}, opts.signals(i))
    end
end

if ~isempty(opts.legendLabels)
    lg = legend(opts.legendLabels);
    lg.Layout.Tile = 'East';
end

postplot(f, opts.filename, 'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true);
end
