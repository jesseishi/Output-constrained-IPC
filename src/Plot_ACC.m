%% Designing a color map.

mycolors = [
     38,  87, 109;
     53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;


%% Part 1A: Step response in laminar wind conditions.
% Reset and load the right things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot\src')

[TTs, Ts_fft, T, S, N, metadata] = loadandprocess2( ...
    ["Results/data/FullDOF_U15_turb0p0_PLExp0p07/*/cIPC_leakyIntegrator_*.csv", ...
     "Results/data/FullDOF_U15_turb0p0_PLExp0p07/refnone/*.csv"]);

model_labels = metadata(2, :);

figDir = 'Results/figures/ACC';


%% Generic check if the simulations look OK.
plotTTsOverview(TTs);


%% Time response in the rotating frame with leaky integrator weight
[f, axs] = preplot(3, 1, 'aspectratio', 1.0667, 'paperFormat', 'ACC', 'column', 1, 'colororder', mycolors);

iLPF = matches('leakyIntegrator', model_labels);
plot(axs(1), TTs{iLPF}(S, :), "FlappingMoment1_MNm")
plot(axs(2), TTs{iLPF}(S, :), "BldPitch1_deg")
plot(axs(3), TTs{iLPF}(S, :), "w_L")
set(axs(3), 'yscale', 'log')

ylabel(axs(1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2), 'pitch (deg)')
ylabel(axs(3), '$$\omega_\mathcal{L}$$', 'interpreter', 'latex')
xlabel('Time (s)')

ylim(axs(2), [11.3, 11.9])
% legend(, 'location', 'northWest')
% legend('\ell^\infty-IPC', '\ell^2-IPC', 'location', 'northWest', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'timeResponseRotatingFrameWithwL.pdf'), 'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1)


%% Time response in the rotating frame withOUT leaky integrator weight
[f, axs] = preplot(2, 1, 'linefrac', 0.33, 'aspectratio', 1.2, 'paperFormat', 'ACC', 'column', 2, 'colororder', mycolors, 'TileSpacing', 'compact');

iLPF = matches('leakyIntegrator', model_labels);
% plot(axs(1), TTs{iLPF}(S, :), "FlappingMoment1_1Pfilt_MNm")
plot(axs(1), TTs{iLPF}(S, :).Time, TTs{iLPF}(S, :).FlappingMoment1_1Pfilt_MNm*1000)
plot(axs(1), TTs{iLPF}(S, :).Time, TTs{iLPF}(S, :).ref_load, '--', 'color', [0.3, 0.3, 0.3])
% plot(axs(1), TTs{iLPF}(S, :), "FlappingMoment1_MNm")
plot(axs(2), TTs{iLPF}(S, :), "BldPitch1_deg")

ylabel(axs(1), {'1P flapping', 'moment (kNm)'})
ylabel(axs(2), 'pitch (deg)')
xlabel('Time (s)')

ylim(axs(1), [-2500, 3500])
yticks(axs(1), -2500:2500:3500)
ylim(axs(2), [11.125, 12.125])
yticks(axs(2), 11.25:0.25:12.25)
xticks(axs(2), seconds(100:40:260))
legend(axs(1), '', 'reference', 'location', 'northeast')
% xticklabels(axs(2), seconds(100:40:260))
xlim(seconds([100, 260]))
postplot(f, fullfile(figDir, 'timeResponseRotatingFrame.pdf'), 'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1);


%% Time response in the nonrotating frame
[f, axs] = preplot(2, 1, 'lineFrac', 0.33, 'aspectratio', 1.2, 'paperFormat', 'ACC', 'column', 2, 'colororder', mycolors);


plot(axs(1), TTs{iLPF}(S, :), "M_ty_mag_filt_kNm")
plot(axs(2), TTs{iLPF}(S, :), "theta_ty_mag_deg")

plot(axs(1), TTs{iLPF}(S, :), "ref_load", 'LineStyle', '--', 'Color', 'k')
ylim(axs(1), [0, 2000])

ylabel(axs(1), {'Tilt/yaw'; 'moment (kNm)'})
ylabel(axs(2), {'Tilt/yaw pitch', 'angle (deg)'})
xlabel('Time (s)')
legend(axs(1), '', 'reference', 'location', 'northeast')
xticks(axs(2), seconds(100:40:260))
xlim(seconds([100, 260]))
postplot(f, fullfile(figDir, 'timeResponseNonrotatingFrame.pdf'), 'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1)


%% Leaky integrator natural frequency
[f, ax] = preplot('lineFrac', 0.33, 'aspectRatio', 1.2, 'paperFormat', 'ACC', 'column', 2, 'YScale', 'log', 'colororder', mycolors);

TTs{iLPF}{S, "w_L"} = max(TTs{iLPF}{S, "w_L"}, 1e-5);
plot(ax, TTs{iLPF}(S, :), "w_L")
ylabel('$$\omega_\mathcal{L}$$', 'interpreter', 'latex')
xlabel('Time (s)')
% xticks(ax, seconds(100:40:260))
% wLmin = 2e-3;
% wLmax = 2e1;
% yline(wLmin, '--', {'min'})
% yline(wLmax, '--', {'max'})
ylim(ax, [1e-3, 1e3])
yticks(logspace(-3, 3, 7))
xticks(seconds(100:40:261))
xlim(seconds([100, 260]))

postplot(f, fullfile(figDir, 'timeResponsewL.pdf'), 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1)

%% Leaky integrator natural frequency and T0
[f, axs] = preplot(2, 1, 'lineFrac', 0.33, 'aspectRatio', 1.2, 'paperFormat', 'ACC', 'column', 2, 'YScale', 'linear', 'colororder', mycolors);

plot(axs(1), TTs{iLPF}(S, :), "T0")

set(axs(2), 'YScale', 'log')
plot(axs(2), TTs{iLPF}(S, :), "w_L")
ylabel(axs(2), '$$\omega_\mathcal{L}$$', 'interpreter', 'latex')
xlabel(axs(2), 'Time (s)')
% xticks(ax, seconds(100:40:260))
% wLmin = 2e-3;
% wLmax = 2e1;
% yline(wLmin, '--', {'min'})
% yline(wLmax, '--', {'max'})
ylim(axs(2), [1e-3, 1e3])
yticks(axs(2), logspace(-3, 3, 3))
xticks(axs(2), seconds(100:40:261))
xlim(axs(2), seconds([100, 260]))

postplot(f, fullfile(figDir, 'timeResponsewLT0.pdf'), 'sharex', true, 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1)

%% Side-by-side
[f, axs] = preplot(2, 3, 'paperFormat', 'ACC', 'column', 2, 'aspectRatio', 2.5, 'colororder', mycolors);
ax_wL = nexttile(3, [2, 1]);

iLPF = matches('leakyIntegrator', model_labels);

plot(axs(1,1), TTs{iLPF}(S, :), "FlappingMoment1_MNm")
% plot(axs(2,1), TTs{iLPF}(S, :), "theta_CPC_deg", 'color', [0.3, 0.3, 0.3]);
% set(axs(2,1), 'colorOrderIndex', 1);
plot(axs(2,1), TTs{iLPF}(S, :), "BldPitch1_deg")
% legend(axs(2,1), 'CPC', 'IPC')
ylabel(axs(1,1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2,1), 'Pitch (deg)')
xlabel(axs(2,1), 'Time (s)')
ylim(axs(2,1), [11.2, 12])
axs(1,1).XTickLabel = [];
axs(1,1).XAxis.Label.String = [];

plot(axs(1,2), TTs{iLPF}(S, :), "M_ty_mag_filt_kNm")
plot(axs(2,2), TTs{iLPF}(S, :), "theta_ty_mag_deg")
plot(axs(1,2), TTs{iLPF}(S, :), "ref_load", 'LineStyle', '--', 'Color', 'k')
ylim(axs(1,2), [0, 2000])
ylabel(axs(1,2), {'Tilt/yaw'; 'moment (kNm)'})
ylabel(axs(2,2), {'Tilt/yaw pitch', 'angle (deg)'})
xlabel(axs(2,2), 'Time (s)')
legend(axs(1,2), '', 'reference', 'location', 'northeast')
axs(1,2).XTickLabel = [];
axs(1,2).XAxis.Label.String = [];

semilogy(ax_wL, TTs{iLPF}(S, :), "w_L")
wLmin = 2e-3;
wLmax = 2e1;
yline(wLmin, '--', {'min'})
yline(wLmax, '--', {'max'})
ylabel(ax_wL, '$\omega_\mathcal{L}$', 'interpreter', 'latex')
grid(ax_wL, 'on');
xlabel('Time (s)')

postplot(f, fullfile(figDir, 'timeResponseCombined.pdf'), 'linkaxes', 'x', 'removeTimetableUnit', true, 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1)

%% On the tilt-yaw plane
[f, axs] = preplot(2, 1, 'linefrac', 0.9, 'aspectratio', 0.8, 'paperFormat', 'ACC', 'column', 1, 'colororder', mycolors);

% Draw no and full IPC (only use the last datapoint)
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end

plot(axs(1), TTs{iLPF}(S, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm")
plot(axs(2), TTs{iLPF}(S, :), "theta_yaw_deg", "theta_tilt_deg")

ylabel(axs(1), 'Tilt moment (kNm)')
xlabel(axs(1), 'Yaw moment (kNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'Yaw pitch (deg)')

% Let's also draw the reference loads in there.
l2refs = unique(TTs{iLPF}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i);
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(axs(1), r*cos(theta), r*sin(theta), 'k--')
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-500, 2500])
ylim(axs(1), [-100, 2000])
xlim(axs(2), [-0.1, 0.5])
ylim(axs(2), [0, 0.4]);

% text(axs(1), -500, 100, 'full IPC')
% text(axs(1), 500, 1350, 'no IPC')
% 
% text(axs(2), -0.12, 0.02, 'no IPC')
% text(axs(2), 0.11, 0.27, 'full IPC')

legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'southeast')
postplot(f, fullfile(figDir, 'ty_response.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L')


%% On the tilt-yaw plane
[f, axs] = preplot(1, 2, 'aspectratio', 1.2, 'paperFormat', 'ACC', 'column', 1);

% Draw no and full IPC (only use the last datapoint)
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end

plot(axs(1), TTs{iLPF}(S, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm")
plot(axs(2), TTs{iLPF}(S, :), "theta_yaw_deg", "theta_tilt_deg")

ylabel(axs(1), 'Tilt moment (kNm)')
xlabel(axs(1), 'yaw moment (kNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'yaw pitch (deg)')

% Let's also draw the reference loads in there.
l2refs = unique(TTs{iLPF}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i);
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(axs(1), r*cos(theta), r*sin(theta), 'k--')
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-500, 1000])
ylim(axs(1), [-100, 2000])
xlim(axs(2), [-0.1, 0.2])
ylim(axs(2), [0, 0.4]);

% text(axs(1), -500, 100, 'full IPC')
% text(axs(1), 500, 1350, 'no IPC')
% 
% text(axs(2), -0.12, 0.02, 'no IPC')
% text(axs(2), 0.11, 0.27, 'full IPC')

lg = legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'southeast');
lg.Layout.Tile = 'North';
postplot(f, fullfile(figDir, 'ty_response.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L')



%% On the tilt-yaw plane
[f, axs] = preplot(1, 2, 'linefrac', 1, 'aspectratio', 1.25, 'paperFormat', 'ACC', 'column', 1);

% Draw no and full IPC (only use the last datapoint)
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end

plot(axs(1), TTs{iLPF}(S, :), "M_yaw_filt_kNm", "M_tilt_filt_kNm")
plot(axs(2), TTs{iLPF}(S, :), "theta_yaw_deg", "theta_tilt_deg")

ylabel(axs(1), 'Tilt moment (kNm)')
xlabel(axs(1), 'Yaw moment (kNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'Yaw pitch (deg)')

% Let's also draw the reference loads in there.
l2refs = unique(TTs{iLPF}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i);
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(axs(1), r*cos(theta), r*sin(theta), 'k--')
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-500, 1000])
ylim(axs(1), [-100, 2000])
xlim(axs(2), [-0.1, 0.2])
ylim(axs(2), [0, 0.4]);

% text(axs(1), -500, 100, 'full IPC')
% text(axs(1), 500, 1350, 'no IPC')
% 
% text(axs(2), -0.12, 0.02, 'no IPC')
% text(axs(2), 0.11, 0.27, 'full IPC')

lg = legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'southeast');
lg.Layout.Tile = 'North';
postplot(f, fullfile(figDir, 'ty_response2.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L')










%% Part 1B: Pareto plot in laminar conditions.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\Matlab utilities'))

% We need all the constant references and no/full IPC.
nameFilter1 = "Results/data/FullDOF_U15_turb8p0_PLExp0p07_2/ref*0/cIPC_adaptive*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
nameFilter2 = "Results/data/FullDOF_U15_turb8p0_PLExp0p07_2/refnone/*.csv";
[TTs, ~, model_labels, ~, ~, T, ~] = loadandprocess([nameFilter1, nameFilter2]);

figDir = 'Results/figures/ACC';

%% Generic check if the simulations look OK.
plotTTsOverview(TTs);

%% Pareto front
[f, axs] = preplot(1, 2, 'column', 2, 'paperFormat', 'ACC', 'aspectRatio', 4);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), T(inoIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), T(inoIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), T(ifullIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), T(ifullIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end

iLPF = contains(model_labels, 'leakyIntegrator');
plotPareto(axs(1), T(iLPF, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean")
plotPareto(axs(2), T(iLPF, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL")

xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (kNm)')
% xlim([0, 1])


legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_laminar.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L')































%% Part 3: turbulent wind conditions.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\Matlab utilities'))

% We need all the constant references and no/full IPC.
% nameFilter1 = "Results/data/FullDOF_U15_TI8_PLexp0p07_seed2*_T2100/*/cIPC_adaptive*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
% nameFilter2 = "Results/data/FullDOF_U15_TI8_PLexp0p07_seed2*_T2100/refnone/IPC*.csv";

namefilter = [...
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed15_T2100/*/*.csv";  % CAUTION! gets all models (ok for now because I've only run it for adaptive and no/full.
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed16_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed17_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed18_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed25_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed26_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed27_T2100/*/*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed28_T2100/*/*.csv"];

tic
[TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(namefilter, 'TStart', seconds(300), 'TEnd', seconds(2100), 'saveTTs', false, 'saveTsfft', false);
% [TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(nameFilter2, 'TStart', seconds(300), 'TEnd', seconds(900), 'saveTTs', true, 'saveTsfft', false);
toc

model_labels = metadata(2, :)';
ref_labels = metadata(3, :)';
seed_labels = metadata(8, :)';

figDir = 'Results/figures/ACC';

% Calculate the average DEL
T.("FlappingMoment_Avg_MNm_DEL") = (T.FlappingMoment1_MNm_DEL + T.FlappingMoment2_MNm_DEL + T.FlappingMoment3_MNm_DEL) / 3;
T.("BldPitch_Avg_deg_DC") = (T.BldPitch1_deg_DC + T.BldPitch2_deg_DC + T.BldPitch3_deg_DC) / 3;

% Add labels
T.("ref_labels") = ref_labels;
T.("model_labels") = model_labels;
T.("seed_labels") = seed_labels;

%% Only run once as backup.
T_full = T;

%% Make a subset of T
idx = contains(seed_labels, ["0", "1", "2", "3", "4", "5"]);
T = T_full(idx, :);
ref_labels = ref_labels(idx);
model_labels = model_labels(idx);
seed_labels = seed_labels(idx);

%% Generic check if the simulations look OK.
% plotTTsOverview(TTs);


%% Pareto front
[f, axs] = preplot(1, 2, 'column', 2, 'paperFormat', 'ACC', 'aspectRatio', 4);

unique_seed_labels = unique(seed_labels);
for iSeed = 1:length(unique_seed_labels)
    seed_label = unique_seed_labels(iSeed);
    seedIdx = matches(seed_labels, seed_label);

    % Add no and full IPC.
    inoIPC = matches(model_labels, 'noIPC') & seedIdx;
    ifullIPC = matches(model_labels, 'fullIPC') & seedIdx;
    scatter(axs(1), T(inoIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'SizeData', 40, 'MarkerEdgeColor', 'k');
    scatter(axs(1), T(ifullIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    scatter(axs(2), T(inoIPC, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'SizeData', 40, 'MarkerEdgeColor', 'k');
    scatter(axs(2), T(ifullIPC, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    for i = 1:length(axs)
        set(axs(i), 'ColorOrderIndex', iSeed)
    end
    
    
    iLPF = contains(model_labels, 'leakyIntegrator') & seedIdx;
    plotPareto(axs(1), T(iLPF, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean")
    plotPareto(axs(2), T(iLPF, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL")
%     break
end

xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (kNm)')
xlim(axs(2), [0, 0.35])


legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');



%% Group and average tryout
% TODO: move this to loadandprocess.
[G, G_ref_labels, G_model_labels] = findgroups(T.ref_labels, T.model_labels);
M_ty_avg = splitapply(@mean, T.M_ty_mag_kNm_mean, G);
theta_ty_avg = splitapply(@mean, T.theta_ty_mag_deg_mean, G);
DEL_avg = splitapply(@mean, T.FlappingMoment_Avg_MNm_DEL, G);
ADC_avg = splitapply(@mean, T.BldPitch_Avg_deg_DC, G);
M_ty_std = splitapply(@std, T.M_ty_mag_kNm_mean, G);  % std is sample standard deviation by default.
theta_ty_std = splitapply(@std, T.theta_ty_mag_deg_mean, G);
DEL_std = splitapply(@std, T.FlappingMoment_Avg_MNm_DEL, G);
ADC_std = splitapply(@std, T.BldPitch_Avg_deg_DC, G);

T_stats = table(theta_ty_avg, M_ty_avg, DEL_avg, ADC_avg, theta_ty_std, M_ty_std, DEL_std, ADC_std);
T_stats

%% Now make a plot of this.
[f, axs] = preplot(1, 2, 'paperFormat', 'ACC', 'column', 2, 'aspectRatio', 4);

GnoIPC = matches(G_model_labels, 'noIPC');
GfullIPC = matches(G_model_labels, 'fullIPC');

scatter(axs(1), T_stats(GnoIPC, :), "theta_ty_avg", "M_ty_avg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), T_stats(GfullIPC, :), "theta_ty_avg", "M_ty_avg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), T_stats(GnoIPC, :), "ADC_avg", "DEL_avg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), T_stats(GfullIPC, :), "ADC_avg", "DEL_avg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end


GLPF = contains(G_model_labels, 'leakyIntegrator');
plotPareto(axs(1), T_stats(GLPF, :), "theta_ty_avg", "M_ty_avg")
plotPareto(axs(2), T_stats(GLPF, :), "ADC_avg", "DEL_avg")

xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (kNm)')
xlim(axs(2), [0, 0.35])


legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent_mean.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');




%% And one with error bars.
[f, axs] = preplot(1, 2, 'paperFormat', 'ACC', 'column', 2, 'aspectRatio', 4);

GnoIPC = matches(G_model_labels, 'noIPC');
GfullIPC = matches(G_model_labels, 'fullIPC');

errorbar(axs(1), T_stats{GnoIPC, "theta_ty_avg"}, T_stats{GnoIPC, "M_ty_avg"}, ...
    T_stats{GnoIPC, "M_ty_std"}, T_stats{GnoIPC, "M_ty_std"}, ...
    T_stats{GnoIPC, "theta_ty_std"}, T_stats{GnoIPC, "theta_ty_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'Color', 'k')
errorbar(axs(1), T_stats{GfullIPC, "theta_ty_avg"}, T_stats{GfullIPC, "M_ty_avg"}, ...
    T_stats{GfullIPC, "M_ty_std"}, T_stats{GfullIPC, "M_ty_std"}, ...
    T_stats{GfullIPC, "theta_ty_std"}, T_stats{GfullIPC, "theta_ty_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')
errorbar(axs(2), T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    T_stats{GnoIPC, "DEL_std"}, T_stats{GnoIPC, "DEL_std"}, ...
    T_stats{GnoIPC, "ADC_std"}, T_stats{GnoIPC, "ADC_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'Color', 'k')
errorbar(axs(2), T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, ...
    T_stats{GfullIPC, "DEL_std"}, T_stats{GfullIPC, "DEL_std"}, ...
    T_stats{GfullIPC, "ADC_std"}, T_stats{GfullIPC, "ADC_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')

for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end


GLPF = contains(G_model_labels, 'leakyIntegrator');
% plotPareto is overrated, it just sorts and then plots. I can do the
% sorting here and then just use a builtin method.
T_stats_sorted = sortrows(T_stats, "theta_ty_avg");
% Here we do a little hack, by setting the sample standard deviation to
% zero on some points we only plot the error bars at certain intervals so
% that the plot doesn't become too messy.
idx = true(length(GLPF), 1);
idx([16, 18]) = false;
T_stats_sorted{GLPF & idx, "theta_ty_std"} = 0;
T_stats_sorted{GLPF & idx, "M_ty_std"} = 0;
T_stats_sorted{GLPF & idx, "ADC_std"} = 0;
T_stats_sorted{GLPF & idx, "DEL_std"} = 0;
errorbar(axs(1), T_stats_sorted{GLPF, "theta_ty_avg"}, T_stats_sorted{GLPF, "M_ty_avg"}, ...
    T_stats_sorted{GLPF, "M_ty_std"}, T_stats_sorted{GLPF, "M_ty_std"}, ...
    T_stats_sorted{GLPF, "theta_ty_std"}, T_stats_sorted{GLPF, "theta_ty_std"})
errorbar(axs(2), T_stats_sorted{GLPF, "ADC_avg"}, T_stats_sorted{GLPF, "DEL_avg"}, ...
    T_stats_sorted{GLPF, "DEL_std"}, T_stats_sorted{GLPF, "DEL_std"}, ...
    T_stats_sorted{GLPF, "ADC_std"}, T_stats_sorted{GLPF, "ADC_std"})

xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (kNm)')
xlim(axs(2), [0, 0.35])

legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent_meanstd.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');




%% Another try with lines.
[f, axs] = preplot(1, 2, 'paperFormat', 'ACC', 'column', 2, 'aspectRatio', 4);

GnoIPC = matches(G_model_labels, 'noIPC');
GfullIPC = matches(G_model_labels, 'fullIPC');

scatter(axs(1), T_stats(GnoIPC, :), "theta_ty_avg", "M_ty_avg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), T_stats(GfullIPC, :), "theta_ty_avg", "M_ty_avg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), T_stats(GnoIPC, :), "ADC_avg", "DEL_avg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), T_stats(GfullIPC, :), "ADC_avg", "DEL_avg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

for i = 1:length(axs)
    set(axs(i), 'ColorOrderIndex', 1)
end

GLPF = contains(G_model_labels, 'leakyIntegrator');
plotPareto(axs(1), T_stats(GLPF, :), "theta_ty_avg", "M_ty_avg")
plotPareto(axs(2), T_stats(GLPF, :), "ADC_avg", "DEL_avg")



[~, isort] = sort(T_stats{:, "theta_ty_avg"});
plot(axs(1), T_stats{isort, "theta_ty_avg"} - T_stats{isort, "theta_ty_std"}, T_stats{isort, "M_ty_avg"} - T_stats{isort, "M_ty_std"})
plot(axs(2), T_stats{isort, "ADC_avg"} - T_stats{isort, "ADC_std"}, T_stats{isort, "DEL_avg"} - T_stats{isort, "DEL_std"})
plot(axs(1), T_stats{isort, "theta_ty_avg"} + T_stats{isort, "theta_ty_std"}, T_stats{isort, "M_ty_avg"} + T_stats{isort, "M_ty_std"})
plot(axs(2), T_stats{isort, "ADC_avg"} + T_stats{isort, "ADC_std"}, T_stats{isort, "DEL_avg"} + T_stats{isort, "DEL_std"})


xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (kNm)')
xlim(axs(2), [0, 0.2])
xlim(axs(1), [-0.01, 0.7])

legend('no IPC', 'full IPC', 'leaky integrator', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent_meanstd.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');

%% Only ADC-DEL.
mycolors = [
     38,  87, 109;
     53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;


[f, ax] = preplot('paperFormat', 'ACC', 'column', 1, 'aspectRatio', 1.75, 'colororder', mycolors);

GnoIPC = matches(G_model_labels, 'noIPC');
GfullIPC = matches(G_model_labels, 'fullIPC');


GLPF = contains(G_model_labels, 'leakyIntegrator');
% plot(ax, T_stats{isort, "ADC_avg"}, T_stats{isort, "DEL_avg"}, 'k')
plot(ax, T_stats{isort, "ADC_avg"}, T_stats{isort, "DEL_avg"}, 'k', 'lineWidth', 1.5)

[~, isort] = sort(T_stats{:, "theta_ty_avg"});
% color = ax.ColorOrder(1, :);
% color = rgb2hsv(color);
% color(3) = 1;%1.5 * color(3);
% color = hsv2rgb(color);
% plot(ax, T_stats{isort, "ADC_avg"} - T_stats{isort, "ADC_std"}, T_stats{isort, "DEL_avg"} - T_stats{isort, "DEL_std"}, 'color', color)
% plot(ax, T_stats{isort, "ADC_avg"} + T_stats{isort, "ADC_std"}, T_stats{isort, "DEL_avg"} + T_stats{isort, "DEL_std"}, 'color', color)

x = min(T_stats{GLPF, "ADC_avg"});
x = [x; T_stats{isort, "ADC_avg"} - T_stats{isort, "ADC_std"}];
x = [x; max(T_stats{GLPF, "ADC_avg"})];
x = [x; T_stats{flipud(isort), "ADC_avg"} + T_stats{flipud(isort), "ADC_std"}];
x = [x; min(T_stats{GLPF, "ADC_avg"})];


y = max(T_stats{GLPF, "DEL_avg"});
y = [y; T_stats{isort, "DEL_avg"} - T_stats{isort, "DEL_std"}];
y = [y; min(T_stats{GLPF, "DEL_avg"})];
y = [y; T_stats{flipud(isort), "DEL_avg"} + T_stats{flipud(isort), "DEL_std"}];
y = [y; max(T_stats{GLPF, "DEL_avg"})];

% patch(ax, x, y, 1, 'FaceAlpha', 0.5, 'FaceColor', mycolors(2,:), 'EdgeColor', [1, 1, 1]);
patchcolor = mycolors(1,:);
patchcolor = rgb2hsv(patchcolor);
patchcolor(3) = 0.7;
patchcolor = hsv2rgb(patchcolor);
patch(ax, x, y, 1, 'FaceColor', patchcolor, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
plot(ax, T_stats{isort, "ADC_avg"}, T_stats{isort, "DEL_avg"}, 'k', 'lineWidth', 1.5)



errorbar(ax, T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    T_stats{GnoIPC, "DEL_std"}, T_stats{GnoIPC, "DEL_std"}, ...
    T_stats{GnoIPC, "ADC_std"}, T_stats{GnoIPC, "ADC_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'Color', 'k')
errorbar(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, ...
    T_stats{GfullIPC, "DEL_std"}, T_stats{GfullIPC, "DEL_std"}, ...
    T_stats{GfullIPC, "ADC_std"}, T_stats{GfullIPC, "ADC_std"}, ...
    'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')


xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, {'Flapping moment damage', 'equivalent load (kNm)'})
xlim(ax, [0, 0.18])
xticks(0:0.02:0.18)

legend('mean', '1-sigma interval', '', 'no IPC', 'full IPC', 'location', 'northeast')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent_meanstd2.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');

