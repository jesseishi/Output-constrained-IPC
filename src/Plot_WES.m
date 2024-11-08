%% Part 1A: Step response in laminar airflow.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% nameFilter = 'Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refstep/*.csv';
[TTs, Ts_fft, T, S, N, metadata] = loadandprocess2( ...
    ["Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refstep/cIPC_l*originalLoad*.csv", ...
    "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/*rnone_Full*.csv"], "TStart", seconds(300), 'TEnd', seconds(425), 'assertEqualLength', false);
% Note that we use a TEnd here. This is because we made the data longer
% than what we're interested in because the 1P filtered load does not look
% good close to the endges of the simulation.

model_labels = metadata(2, :);
figDir = 'Results/figures/WES';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

% Change the order of the timetables, so that linf comes first.
newIdx = [2, 1, 3, 4];
TTs = TTs(newIdx);
model_labels = model_labels(newIdx);

originalLoadIndices = find(contains(model_labels, 'originalLoad'))';

%% Generic check if the simulations look OK.
plotTTsOverview(TTs, 'legendLabels', model_labels)

%% Time response in the rotating frame
[f, axs] = preplot(3, 1, 'paperFormat', 'WES', 'column', 1, 'colororder', mycolors, 'aspectRatio', 1);

for j = 1:length(originalLoadIndices)
    i = originalLoadIndices(j);
    plot(axs(1), TTs{i}(S, :), "FlappingMoment1_MNm")
    
    plot(axs(2), TTs{i}(S, :), "FlappingMoment1_1Pfilt_MNm")
    if j == length(originalLoadIndices)
        plot(axs(2), TTs{i}(S, :).Time, TTs{i}(S, :).ref_load ./ 1e3, 'linestyle', '--', 'color', [0.3, 0.3, 0.3])
    end
    plot(axs(3), TTs{i}(S, :), "BldPitch1_deg")
    %     plot(axs(2), TTs{i}(S, :), "theta_1", 'linestyle', '--')  % Looks exactly the same = good.
end

ylim(axs(1), [28, 34])
yticks(axs(1), 26:2:34)

ylim(axs(2), [-2, 2])
yticks(axs(2), -2:1:2)

% % I don't like this but it works.
% lines = findobj(axs(2), 'Type', 'Line');
% legend(axs(1), lines(end:-1:1), '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');

lg = legend(axs(2), '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'interpreter', 'latex');
lg.Layout.Tile = 'south';
% lg.Orientation = 'horizontal';

ylim(axs(3), [11.25, 12.0])
yticks(axs(3), 11.0:0.25:12.0)

xticks(axs(1), seconds(300:25:425))
xticks(axs(2), seconds(300:25:425))
xticks(axs(3), seconds(300:25:425))

ylabel(axs(1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2), {'1P flapping', 'moment (MNm)'})
ylabel(axs(3), 'Pitch angle (deg)')
xlabel(axs(3), 'Time (s)')
% legend(model_labels(originalLoadIndices), 'location', 'northWest')
% legend('\ell^\infty-IPC', '\ell^2-IPC', 'location', 'northWest', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'timeResponseRotatingFrame.pdf'), ...
    'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1);



%% Time response in the nonrotating frame
[f, axs] = preplot(2, 1, 'paperFormat', 'WES', 'column', 1, 'colororder', mycolors);

for j = 1:length(originalLoadIndices)
    i = originalLoadIndices(j);
    plot(axs(1), TTs{i}(S, :).Time, TTs{i}(S, :).M_ty_mag_filt_kNm/1e3)
    plot(axs(2), TTs{i}(S, :), "theta_ty_mag_deg")
end

i_l2 = find(contains(model_labels, 'l2'));
plot(axs(1), TTs{i_l2}(S, :).Time, TTs{i_l2}(S, :).ref_load/1e3, 'LineStyle', '--', 'Color', 'k')
ylim(axs(1), [0, 2])

ylabel(axs(1), {'nonrotating', 'moment (MNm)'})
ylabel(axs(2), {'nonrotating', 'pitch angle (deg)'})
xlabel('Time (s)')

xticks(axs(2), seconds(300:25:425))
xticks(axs(1), seconds(300:25:425))

% I don't like this but it works.
lines = findobj(axs(1), 'Type', 'Line');
legend(axs(2), lines(end:-1:1), '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');
% legend(axs(1), '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'northeast', 'interpreter', 'latex')

postplot(f, fullfile(figDir, 'timeResponseNonrotatingFrame.pdf'), 'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1);


%% On the tilt-yaw plane
[f, axs] = preplot(2, 1, 'paperFormat', 'WES', 'column', 1, 'colororder', mycolors, 'aspectRatio', 0.73);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

for i = 1:2
    set(axs(i), 'colororderindex', 1);
end

for j = 1:length(originalLoadIndices)
    i = originalLoadIndices(j);
    plot(axs(1), TTs{i}(S, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm")
    plot(axs(2), TTs{i}(S, :), "theta_yaw_deg", "theta_tilt_deg")
end

ylabel(axs(1), 'Tilt moment (MNm)')
xlabel(axs(1), 'yaw moment (MNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'yaw pitch (deg)')

% Let's also draw the reference loads in there.
color = mycolors(2, :);
color = rgb2hsv(color);
color(2) = 0.5 * color(2);  % reduce chroma.
color = hsv2rgb(color);

l2refs = unique(TTs{find(contains(model_labels, "l2"), 1)}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i) ./ 1e3;
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(axs(1), r*cos(theta), r*sin(theta), '--', 'color', color)
end

color = mycolors(1, :);
color = rgb2hsv(color);
color(2) = 0.5 * color(2);  % reduce chroma.
color = hsv2rgb(color);

linfrefs = unique(TTs{find(contains(model_labels, "linf"), 1)}.ref_load);
for i = 1:length(linfrefs)
    r = linfrefs(i) ./1e3 ;
    % Draw a box.
    x = [r, -r, -r, r, r];
    y = [r, r, -r, -r, r];
    plot(axs(1), x, y, '--', 'color', color)
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-1, 2])
ylim(axs(1), [-0.1, 2])
xlim(axs(2), [-0.2, 0.4])
ylim(axs(2), [0, 0.4]);

% legend(['no IPC'; 'full IPC'; model_labels(originalLoadIndices)], 'location', 'southeast')
% lg = legend(axs(1), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'northeast', 'interpreter', 'latex');
% lg.Layout.Tile = 'East';

% This is also not idea.
plot(axs(2), -10:-9, 0:1, '--k')  % Phantom line.
legend(axs(2), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'ty_response.pdf'), ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1)

























%% Part 1B: Pareto plot in laminar conditions.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% We need all the constant references and no/full IPC.
nameFilter1 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/ref*0/cIPC_l*originalLoad*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
nameFilter2 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/*rnone_Full*.csv";
[~, ~, T, ~, ~, metadata] = loadandprocess2([nameFilter1, nameFilter2], ...
    'assertEqualLength', false, 'TStart', seconds(300), 'TEnd', seconds(400), 'saveTTs', false, 'saveTsfft', false);
% Not all simulations were the same length -> Only use the part from 300 to
% 400 seconds. This is long enough since we're doing laminar flow.

figDir = 'Results/figures/WES';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

model_labels = metadata(2, :);

fprintf('Done.\n');

%% Generic check if the simulations look OK.
% plotTTsOverview(TTs)

%% Pareto front
[f, axs] = preplot(1, 2, 'lineFrac', 1.0, 'paperFormat', 'WES', 'aspectRatio', 2.5, 'colororder', mycolors, 'column', 2);

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


unique_model_labels = unique(model_labels);
unique_cIPC_model_labels = unique_model_labels(~contains(unique_model_labels, 'IPC'));
unique_cIPC_model_labels = unique_cIPC_model_labels([2, 1]);  % Change the order so that linf comes first.
for i = 1:length(unique_cIPC_model_labels)
    model_label = unique_cIPC_model_labels(i);
    idx = matches(model_labels, model_label);
    
    plotPareto(axs(1), T(idx, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'linestyle', '-x')
    %     plotPareto(axs(2), T(idx, :), "BldPitch1_deg_psd_at1P", "FlappingMoment1_MNm_psd_at1P")
    plotPareto(axs(2), T(idx, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'linestyle', '-x')
    %     disp(idx)
end

% NOTE: there's quite a sharp transition for Linf, not that this is
% supposed to be there because we transition from only tilt to tilt and yaw
% at 45 degrees.


xlabel(axs(1), 'Mean nonrotating pitch (deg)')
ylabel(axs(1), 'Mean nonrotating moment (MNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment damage equivalent load (MNm)')
% xlim([0, 1])


legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'tradeoff_laminar.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1)

%% Only DEL vs ADC
[f, ax] = preplot(1, 1, 'lineFrac', 1.0, 'paperFormat', 'WES', 'aspectRatio', 1.618, 'colororder', mycolors, 'column', 1);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, T(inoIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(ax, T(ifullIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

set(ax, 'ColorOrderIndex', 1)


unique_model_labels = unique(model_labels);
unique_cIPC_model_labels = unique_model_labels(~contains(unique_model_labels, 'IPC'));
unique_cIPC_model_labels = unique_cIPC_model_labels([2, 1]);  % Change the order so that linf comes first.
for i = 1:length(unique_cIPC_model_labels)
    model_label = unique_cIPC_model_labels(i);
    idx = matches(model_labels, model_label);
    
    plotPareto(ax, T(idx, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL")
    %     disp(idx)
end



xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, {'Flapping moment', 'damage equivalent load (MNm)'})

legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_laminar.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1)


%% Also calculate percentage-wise the gains.
iInf = matches(model_labels, 'linforiginalLoad');
iL2 = matches(model_labels, 'l2originalLoad');

xLInf = T{iInf, "theta_ty_mag_deg_mean"};
yLInf = T{iInf, "M_ty_mag_kNm_mean"};
xL2 = T{iL2, "theta_ty_mag_deg_mean"};
yL2 = T{iL2, "M_ty_mag_kNm_mean"};

x = linspace(0, max([xLInf; xL2]), 1e3);
method = 'linear';
y1 = interp1(xLInf, yLInf, x, method);
y2 = interp1(xL2, yL2, x, method);

y2wrty1 = (y2 - y1) ./ y1;

preplot;
plot(x, y2wrty1)

xLInf = T{iInf, "BldPitch1_deg_DC"};
yLInf = T{iInf, "FlappingMoment1_MNm_DEL"};
xL2 = T{iL2, "BldPitch1_deg_DC"};
yL2 = T{iL2, "FlappingMoment1_MNm_DEL"};

x = linspace(0, max([xLInf; xL2]), 1e3);
y1 = interp1(xLInf, yLInf, x, method);  % TODO: also try spline and pchip interpolation though there should be no difference once we get enough datapoints.
y2 = interp1(xL2, yL2, x, method);

y2wrty1 = (y2 - y1) ./ y1;
plot(x, y2wrty1)

legend('mean moment vs mean pitch', 'DEL vs ADC')
title('L2 w.r.t. Linf')

postplot;














%% Part 1B: Pareto plot in laminar conditions with a distribution for full IPC.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% We need all the constant references and no/full IPC.
nameFilter1 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/ref*0/cIPC_l*originalLoad*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
nameFilter2 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/IPC_fullIPC_rnone_FullDOF*.csv";
[~, ~, T, ~, ~, metadata] = loadandprocess2([nameFilter1, nameFilter2], ...
    'assertEqualLength', false, 'TStart', seconds(300), 'TEnd', seconds(400), 'saveTTs', false, 'saveTsfft', false);
% Not all simulations were the same length -> Only use the part from 300 to
% 400 seconds. This is long enough since we're doing laminar flow.
nameFilter = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/IPC_fullIPC_rnone_w*.csv";
[~, ~, T2, ~, ~, metadata2] = loadandprocess2(nameFilter, ...
    'assertEqualLength', false, 'TStart', seconds(300), 'TEnd', seconds(400), 'saveTTs', false, 'saveTsfft', false);


figDir = 'Results/figures/WES';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

model_labels = metadata(2, :);

fprintf('Done.\n');


%% Calculate average over three blades.
T2.("FlappingMoment_Avg_MNm_DEL") = (T2.FlappingMoment1_MNm_DEL + T2.FlappingMoment2_MNm_DEL + T2.FlappingMoment3_MNm_DEL) / 3;
T2.("BldPitch_Avg_deg_DC") = (T2.BldPitch1_deg_DC + T2.BldPitch2_deg_DC + T2.BldPitch3_deg_DC) / 3;

% So the variance is really too small here to put in a plot.
mean(T2.FlappingMoment_Avg_MNm_DEL)
std(T2.FlappingMoment_Avg_MNm_DEL)
mean(T2.BldPitch_Avg_deg_DC)
std(T2.BldPitch_Avg_deg_DC)

%% Generic check if the simulations look OK.
% plotTTsOverview(TTs)

%% Pareto front
[f, axs] = preplot(1, 2, 'lineFrac', 1.0, 'paperFormat', 'WES', 'aspectRatio', 2.5, 'colororder', mycolors, 'column', 2);

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


unique_model_labels = unique(model_labels);
unique_cIPC_model_labels = unique_model_labels(~contains(unique_model_labels, 'IPC'));
unique_cIPC_model_labels = unique_cIPC_model_labels([2, 1]);  % Change the order so that linf comes first.
for i = 1:length(unique_cIPC_model_labels)
    model_label = unique_cIPC_model_labels(i);
    idx = matches(model_labels, model_label);
    
    plotPareto(axs(1), T(idx, :), "theta_ty_mag_deg_mean", "M_ty_mag_kNm_mean", 'linestyle', '-x')
    %     plotPareto(axs(2), T(idx, :), "BldPitch1_deg_psd_at1P", "FlappingMoment1_MNm_psd_at1P")
    plotPareto(axs(2), T(idx, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'linestyle', '-x')
    %     disp(idx)
end

% NOTE: there's quite a sharp transition for Linf, not that this is
% supposed to be there because we transition from only tilt to tilt and yaw
% at 45 degrees.


xlabel(axs(1), 'Mean nonrotating pitch (deg)')
ylabel(axs(1), 'Mean nonrotating moment (kNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((kNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment damage equivalent load (kNm)')
% xlim([0, 1])


legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'tradeoff_laminar.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1)

%% Only DEL vs ADC
[f, ax] = preplot(1, 1, 'lineFrac', 1.0, 'paperFormat', 'WES', 'aspectRatio', 1.618, 'colororder', mycolors, 'column', 1);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, T(inoIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'SizeData', 40, 'MarkerEdgeColor', 'k');
scatter(ax, T(ifullIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

set(ax, 'ColorOrderIndex', 1)


unique_model_labels = unique(model_labels);
unique_cIPC_model_labels = unique_model_labels(~contains(unique_model_labels, 'IPC'));
unique_cIPC_model_labels = unique_cIPC_model_labels([2, 1]);  % Change the order so that linf comes first.
for i = 1:length(unique_cIPC_model_labels)
    model_label = unique_cIPC_model_labels(i);
    idx = matches(model_labels, model_label);
    
    plotPareto(ax, T(idx, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL")
    %     disp(idx)
end



xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, {'Flapping moment', 'damage equivalent load (kNm)'})

legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_laminar.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1)


%% Also calculate percentage-wise the gains.
iInf = matches(model_labels, 'linforiginalLoad');
iL2 = matches(model_labels, 'l2originalLoad');

xLInf = T{iInf, "theta_ty_mag_deg_mean"};
yLInf = T{iInf, "M_ty_mag_kNm_mean"};
xL2 = T{iL2, "theta_ty_mag_deg_mean"};
yL2 = T{iL2, "M_ty_mag_kNm_mean"};

x = linspace(0, max([xLInf; xL2]), 1e3);
method = 'linear';
y1 = interp1(xLInf, yLInf, x, method);
y2 = interp1(xL2, yL2, x, method);

y2wrty1 = (y2 - y1) ./ y1;

preplot;
plot(x, y2wrty1)

xLInf = T{iInf, "BldPitch1_deg_DC"};
yLInf = T{iInf, "FlappingMoment1_MNm_DEL"};
xL2 = T{iL2, "BldPitch1_deg_DC"};
yL2 = T{iL2, "FlappingMoment1_MNm_DEL"};

x = linspace(0, max([xLInf; xL2]), 1e3);
y1 = interp1(xLInf, yLInf, x, method);  % TODO: also try spline and pchip interpolation though there should be no difference once we get enough datapoints.
y2 = interp1(xL2, yL2, x, method);

y2wrty1 = (y2 - y1) ./ y1;
plot(x, y2wrty1)

legend('mean moment vs mean pitch', 'DEL vs ADC')
title('L2 w.r.t. Linf')

postplot;




















%% Part 2: rotating wind but with a different reference.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% We need all the constant references and no/full IPC.
nameFilter1 = "Results/data/FullDOF_U15_RotatingWind_PLExp0p07_T900/ref1000/cIPC_l*originalLoad*_w0.2_wo2.0_*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
nameFilter2 = "Results/data/FullDOF_U15_RotatingWind_PLExp0p07_T900/refnone/*_w0.2_wo2.0_*.csv";
[TTs, ~, ~, S, N, metadata] = loadandprocess2([nameFilter1, nameFilter2], 'TStart', seconds(300), 'TEnd', seconds(420));

model_labels = metadata(2, :);
ref_labels = metadata(3, :);
figDir = 'Results/figures/WES';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

% Change the order of the timetables, so that linf comes first.
newIdx = [2, 1, 3, 4];
TTs = TTs(newIdx);
model_labels = model_labels(newIdx);
ref_labels = ref_labels(newIdx);


idxNoIPC = contains(model_labels, 'noIPC');
idxFullIPC = contains(model_labels, 'fullIPC');
idxLinf = contains(model_labels, 'linf');
idxL2 = contains(model_labels, 'l2');

%% General check
plotTTsOverview(TTs, "legendLabels", model_labels)
plotTTsOverview(TTs, "legendLabels", model_labels, signals=["M_tilt_MNm", "M_yaw_MNm"])


%% For the Linf controller integrator saturation and reference load sign are important.
% Need to show: tilt and yaw load and tilt and yaw pitch.
[f, axs] = preplot(3, 1, 'paperFormat', 'WES', 'aspectRatio', 0.95, 'colororder', mycolors, 'column', 1);

plot(axs(1), TTs{idxLinf}(S, :).Time, TTs{idxLinf}(S, :).ref_load_tilt./1e3, 'color', 'black', 'linestyle', '--')
set(axs(1), 'ColorOrderIndex', 1)
plot(axs(1), TTs{idxLinf}(S, :), "M_tilt_filt_MNm")
ylabel(axs(1), 'Tilt moment (MNm)')
legend(axs(1), 'Reference', '$\ell^\infty$-IPC', 'location', 'southwest', 'interpreter', 'latex')
ylim(axs(1), [-2, 2])
yticks(axs(1), -2:1:2)

plot(axs(2), TTs{idxNoIPC}(S, :), "M_tilt_filt_MNm", 'color', 'black', 'linestyle', '--')
set(axs(2), 'ColorOrderIndex', 1)
plot(axs(2), TTs{idxLinf}(S, :), "M_tilt_0_MNm")
ylabel(axs(2), {'Original', 'tilt moment (MNm)'})
legend(axs(2), 'No IPC', '$\ell^\infty$-IPC', 'interpreter', 'latex', 'location', 'southwest')
ylim(axs(2), [-2, 2])
yticks(axs(2), -2:1:2)

plot(axs(3), TTs{idxLinf}(S, :), "theta_tilt_deg")
ylabel(axs(3), 'Tilt pitch (deg)')
ylim(axs(3), [-0.2, 0.2])
yticks(axs(3), -0.2:0.1:0.2)
xlabel(axs(3), 'Time (s)')

postplot(f, fullfile(figDir, 'rotatingWindLinf.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1, 'removeTimetableUnit', true, 'sharex', true);

%% For the L2 controller need to show the phase of the pitch signal (should be ideal)
[f, ax] = preplot(1, 1, 'paperFormat', 'WES', 'aspectRatio', 2, 'colororder', mycolors, 'column', 1);

% plot(axs(1), TTs{idxL2}(S, :), "M_ty_0_mag_kNm")
% plot(axs(1), TTs{idxNoIPC}(S, :), "M_ty_mag_filt_kNm", 'color', 'black', 'linestyle', '--')

% plot(axs(2), TTs{idxL2}(S, :), "M_ty_phase_deg")
% plot(axs(2), TTs{idxL2}(S, :), "M_ty_0_phase_deg")
% plot(axs(2), TTs{idxFullIPC}(S, :), "M_ty_phase_deg", 'color', 'black', 'linestyle', '--')

set(ax, 'ColorOrderIndex', 2)
plot(ax, TTs{idxL2}(S, :), "M_ty_0_phase_deg")  % Not the pitch signal phase signal but that one is difficult to filter, but theta_phase == M_phase for the l2 controller.
plot(ax, TTs{idxFullIPC}(S, :), "theta_ty_phase_deg", 'color', 'black', 'linestyle', '--')
legend('$\ell^2$-IPC', 'full IPC', 'interpreter', 'latex')
ylabel(ax, 'Pitch signal phase (deg)')
xlabel(ax, 'Time (s)')

ylim(ax, [-90, 90])
yticks(ax, -90:45:90)

postplot(f, fullfile(figDir, 'rotatingWindL2.pdf'), 'fontname', 'NimbusRomNo9L', 'fontsize', 8, 'linewidth', 1, 'removeTimetableUnit', true)























%% Part 3: turbulent wind conditions.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

namefilter = [...
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed100_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed101_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed102_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed103_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed104_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed105_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed106_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed107_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed108_T2100/*/*_wo2.0*.csv";
    "Results/data/FullDOF_U15_TI8_PLExp0p07_seed109_T2100/*/*_wo2.0*.csv"];
tic
[~, ~, T, S, N, metadata] = loadandprocess2(namefilter, 'TStart', seconds(300), 'TEnd', seconds(2100), 'saveTTs', false, 'saveTsfft', false);
% [TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(nameFilter2, 'TStart', seconds(300), 'TEnd', seconds(900), 'saveTTs', true, 'saveTsfft', false);
toc

model_labels = metadata(2, :)';
ref_labels = metadata(3, :)';
seed_labels = metadata(10, :)';
gain_labels = metadata(4, :)';


figDir = 'Results/figures/WES';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

% TODO: move to loadandprocess2
% Calculate the average DEL and ADC
T.("FlappingMoment_Avg_MNm_DEL") = (T.FlappingMoment1_MNm_DEL + T.FlappingMoment2_MNm_DEL + T.FlappingMoment3_MNm_DEL) / 3;
T.("BldPitch_Avg_deg_DC") = (T.BldPitch1_deg_DC + T.BldPitch2_deg_DC + T.BldPitch3_deg_DC) / 3;

% Add labels
T.("ref_labels") = ref_labels;
T.("model_labels") = model_labels;
T.("seed_labels") = seed_labels;
T.("gain_labels") = gain_labels;

% Once as backup.
save('src/TempCache/T.mat', 'T')
T_full = T;

beep

%% Overview of all the different seeds
[f, axs] = preplot(1, 2, 'column', 2, 'paperFormat', 'WES', 'aspectRatio', 4);

unique_seed_labels = unique(seed_labels);
for iSeed = 1:length(unique_seed_labels)
    seed_label = unique_seed_labels(iSeed);
    seedIdx = matches(seed_labels, seed_label);
    
    % Add no and full IPC.
    inoIPC = matches(model_labels, 'noIPC') & seedIdx;
    ifullIPC = matches(model_labels, 'fullIPC') & seedIdx;
    scatter(axs(1), T(inoIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'SizeData', 40, 'MarkerEdgeColor', 'k');
    scatter(axs(1), T(ifullIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    scatter(axs(2), T(inoIPC, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'SizeData', 40, 'MarkerEdgeColor', 'k');
    scatter(axs(2), T(ifullIPC, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    for i = 1:length(axs)
        set(axs(i), 'ColorOrderIndex', iSeed)
    end
    
    
    iLPF = contains(model_labels, 'l2') & seedIdx;
    plotPareto(axs(1), T(iLPF, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'linestyle', 'x-')
    plotPareto(axs(2), T(iLPF, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'linestyle', 'x-')
    for i = 1:length(axs)
        set(axs(i), 'ColorOrderIndex', iSeed)
    end
    iLPF = contains(model_labels, 'linf') & seedIdx;
    plotPareto(axs(1), T(iLPF, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'linestyle', 'x--')
    plotPareto(axs(2), T(iLPF, :), "BldPitch_Avg_deg_DC", "FlappingMoment_Avg_MNm_DEL", 'linestyle', 'x--')
    %     break
end

xlabel(axs(1), 'Mean tilt-yaw pitch (deg)')
ylabel(axs(1), 'Mean tilt/yaw moment (MNm)')
% xlabel(axs(2), 'Pitch 1P power (deg$^2$)')
% ylabel(axs(2), 'Flapping moment 1P power ((MNm)$^2$)')
xlabel(axs(2), 'Actuator duty cycle (-)')
ylabel(axs(2), 'Flapping moment DEL (MNm)')
xlim(axs(2), [0, 0.22])


legend('no IPC', 'full IPC', '$\ell^2$-IPC', '$\ell^\infty$-IPC', 'interpreter', 'latex', 'location', 'northeast')
postplot(f, fullfile(figDir, 'tradeoffTurbulentOverview.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');



%% Group and average tryout
% TODO: move this to loadandprocess.
[G, G_ref_labels, G_model_labels, G_gain_labels] = findgroups(T.ref_labels, T.model_labels, T.gain_labels);
M_ty_avg = splitapply(@mean, T.M_ty_mag_MNm_mean, G);
theta_ty_avg = splitapply(@mean, T.theta_ty_mag_deg_mean, G);
DEL_avg = splitapply(@mean, T.FlappingMoment_Avg_MNm_DEL, G);
ADC_avg = splitapply(@mean, T.BldPitch_Avg_deg_DC, G);
M_ty_std = splitapply(@std, T.M_ty_mag_MNm_mean, G);  % std is sample standard deviation by default.
theta_ty_std = splitapply(@std, T.theta_ty_mag_deg_mean, G);
DEL_std = splitapply(@std, T.FlappingMoment_Avg_MNm_DEL, G);
ADC_std = splitapply(@std, T.BldPitch_Avg_deg_DC, G);

T_stats = table(theta_ty_avg, M_ty_avg, DEL_avg, ADC_avg, theta_ty_std, M_ty_std, DEL_std, ADC_std);
T_stats.ref_labels = G_ref_labels;
T_stats.model_labels = G_model_labels;
T_stats.gain_labels = G_gain_labels;

% Remove certain gains.
iSelect = ~matches(T_stats.gain_labels, "w1.50");
T_stats = T_stats(iSelect, :);

[~, iSort] = sort(T_stats.ADC_avg);
T_stats = T_stats(iSort, :);
T_stats


%% DEL vs ADC for conventional IPC.
[f, ax] = preplot(1, 1, 'fnum', 2, 'aspectRatio', 1.25, 'column', 1, 'paperFormat', 'WES');

GnoIPC = matches(T_stats.model_labels, 'noIPC');
GfullIPC = matches(T_stats.model_labels, 'fullIPC');


% Full IPC patch
alpha = 0.15;
drawStandardDeviationPatch(ax, T_stats, GfullIPC, [0, 0, 0], alpha, 0);
plot(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, '-xk')

% Full IPC chosen in this work
idx = matches(T_stats.gain_labels, 'w0.2') & matches(T_stats.model_labels, 'fullIPC');
% scatter(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, 75, 'k', 'linewidth', 1)
% errorbar(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, ...
%     T_stats{idx, "DEL_std"}, T_stats{idx, "DEL_std"}, ...
%     T_stats{idx, "ADC_std"}, T_stats{idx, "ADC_std"}, ...
%     'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'color', 'k')

% no IPC
markersize = 5;  % 10 for ppt.
errorbar(ax, T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    T_stats{GnoIPC, "DEL_std"}, T_stats{GnoIPC, "DEL_std"}, ...
    T_stats{GnoIPC, "ADC_std"}, T_stats{GnoIPC, "ADC_std"}, ...
    'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'Color', 'k')

% Labels and stuff
xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, 'Damage equivalent load (kNm)')
xlim(ax, [0, 0.2])
xticks(0:0.05:0.2)
ylim([11, 15]);

% legend('', 'full IPC', '$w_c = 0.2$ rad/s', 'no IPC', 'location', 'northeast', 'interpreter', 'latex')
legend('', 'full IPC', 'no IPC', 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent_conventional.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1);


%% DEL vs ADC for L2 and Linf IPC.
[f, ax] = preplot('paperFormat', 'WES', 'column', 1, 'aspectRatio', 1.25, 'colororder', mycolors);
% [f, ax] = preplot('paperFormat', 'ppt', 'aspectRatio', 1.618, 'linefrac', 0.65, 'colororder', mycolors);
markersize = 5;  % 10 for ppt.


GLinf = matches(T_stats.model_labels, 'linforiginalLoad');
GL2 = matches(T_stats.model_labels, 'l2originalLoad');
GnoIPC = matches(T_stats.model_labels, 'noIPC');
GfullIPC = matches(T_stats.model_labels, 'fullIPC');

% Linf patch
drawStandardDeviationPatch(ax, T_stats, GLinf, mycolors(1, :), alpha, 0.7);
plot(ax, T_stats{GLinf, "ADC_avg"}, T_stats{GLinf, "DEL_avg"}, '-x', 'color', mycolors(1,:))


% L2 patch
drawStandardDeviationPatch(ax, T_stats, GL2, mycolors(3, :), alpha, 0.7);
plot(ax, T_stats{GL2, "ADC_avg"}, T_stats{GL2, "DEL_avg"}, '-x', 'color', mycolors(2,:))



% errorbar(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, ...
%     T_stats{GfullIPC, "DEL_std"}, T_stats{GfullIPC, "DEL_std"}, ...
%     T_stats{GfullIPC, "ADC_std"}, T_stats{GfullIPC, "ADC_std"}, ...
%     'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')
% Full IPC patch
drawStandardDeviationPatch(ax, T_stats, GfullIPC, [0, 0, 0], alpha, 0);
plot(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, '-xk')
scatter(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, 75, 'filled', 'k', 'linewidth', 1)
% errorbar(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, ...
%     T_stats{idx, "DEL_std"}, T_stats{idx, "DEL_std"}, ...
%     T_stats{idx, "ADC_std"}, T_stats{idx, "ADC_std"}, ...
%     'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'color', 'k')
errorbar(ax, T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    T_stats{GnoIPC, "DEL_std"}, T_stats{GnoIPC, "DEL_std"}, ...
    T_stats{GnoIPC, "ADC_std"}, T_stats{GnoIPC, "ADC_std"}, ...
    'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'Color', 'k')


xlim(ax, [0, 0.2])
xticks(0:0.05:0.2)
xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, 'Damage equivalent load (MNm)')

legend('', '$\ell^\infty$-IPC', '', '$\ell^2$-IPC', '', 'full IPC', '$w_c = 0.2$ rad/s', 'no IPC', 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_turbulent.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1);

% pptDir = replace(figDir, 'WES', 'NAWEA');
% postplot(f, fullfile(pptDir, 'DELvsADC_turbulent.svg'), 'fontsize', 15, 'fontname', 'NimbusRomNo9L');

%% DEL vs ADC for L2 and Linf IPC for presentation
% [f, ax] = preplot('paperFormat', 'WES', 'column', 1, 'aspectRatio', 1.25, 'colororder', mycolors);
[f, ax] = preplot('paperFormat', 'ppt', 'aspectRatio', 1.618, 'linefrac', 0.65, 'colororder', mycolors);
markersize = 10;  % 10 for ppt.


GLinf = matches(T_stats.model_labels, 'linforiginalLoad');
GL2 = matches(T_stats.model_labels, 'l2originalLoad');
GnoIPC = matches(T_stats.model_labels, 'noIPC');
GfullIPC = matches(T_stats.model_labels, 'fullIPC');

alpha = 0.15;

% Linf patch
drawStandardDeviationPatch(ax, T_stats, GLinf, mycolors(1, :), alpha, 0.7);
plot(ax, T_stats{GLinf, "ADC_avg"}, T_stats{GLinf, "DEL_avg"}, '-x', 'color', mycolors(1,:), 'lineWidth', 1.5)


% L2 patch
drawStandardDeviationPatch(ax, T_stats, GL2, mycolors(3, :), alpha, 0.7);
plot(ax, T_stats{GL2, "ADC_avg"}, T_stats{GL2, "DEL_avg"}, '-x', 'color', mycolors(2,:), 'lineWidth', 1.5)


errorbar(ax, T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    T_stats{GnoIPC, "DEL_std"}, T_stats{GnoIPC, "DEL_std"}, ...
    T_stats{GnoIPC, "ADC_std"}, T_stats{GnoIPC, "ADC_std"}, ...
    'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'Color', 'k')
% errorbar(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, ...
%     T_stats{idx, "DEL_std"}, T_stats{idx, "DEL_std"}, ...
%     T_stats{idx, "ADC_std"}, T_stats{idx, "ADC_std"}, ...
%     'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')
% Full IPC patch
% drawStandardDeviationPatch(ax, T_stats, GfullIPC, [0, 0, 0], alpha, 0);
% plot(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, '-xk', 'lineWidth', 1.5)
errorbar(ax, T_stats{idx, "ADC_avg"}, T_stats{idx, "DEL_avg"}, ...
    T_stats{idx, "DEL_std"}, T_stats{idx, "DEL_std"}, ...
    T_stats{idx, "ADC_std"}, T_stats{idx, "ADC_std"}, ...
    'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'color', 'k')

xlim(ax, [0, 0.18])
xticks(0:0.05:0.22)
xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, 'Damage equivalent load (MNm)')

legend('', '$\ell^\infty$-IPC', '', '$\ell^2$-IPC', 'no IPC', '', 'full IPC', 'location', 'northeast', 'interpreter', 'latex')
% postplot(f, fullfile(figDir, 'DELvsADC_turbulent.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L');

pptDir = replace(figDir, 'WES', 'NAWEA');
postplot(f, fullfile(pptDir, 'DELvsADC_turbulent.svg'), 'fontsize', 15, 'fontname', 'NimbusRomNo9L');


%% 'Lorenz' curve (called P-P plot, or percent-percent plot) - I think: need to check.
% I think this should be in there, maybe also with standard deviations.
% For 5 (!) turbulent wind conditions, it looks pretty crazy, we can get
% 90% of the DEL reduction with 50% of the ADC increase.
[f, ax] = preplot('paperFormat', 'WES', 'column', 1, 'aspectRatio', 1, 'linefrac', 0.8, 'colororder', mycolors);
% [f, ax] = preplot('paperFormat', 'ppt', 'linefrac', 0.65, 'aspectRatio', 1.618, 'colororder', mycolors);
markersize = 5;

Gs = [GLinf, GL2, GfullIPC];
% xs = [xLinf, xL2];
% ys = [yLinf, yL2];
patchcolors = [mycolors(1,:); mycolors(2,:); [0,0,0]];
transparencies = [0.7, 0.7, 0];

T_stats_normalized = T_stats;
T_stats_normalized.ADC_avg = 100*(T_stats.ADC_avg - T_stats{GnoIPC, "ADC_avg"}) ./ (T_stats{idx, "ADC_avg"} - T_stats{GnoIPC, "ADC_avg"});
T_stats_normalized.DEL_avg = 100*(T_stats.DEL_avg - T_stats{GnoIPC, "DEL_avg"}) ./ (T_stats{idx, "DEL_avg"} - T_stats{GnoIPC, "DEL_avg"});

T_stats_normalized.ADC_std = 100*(T_stats.ADC_std - 0*T_stats{GnoIPC, "ADC_std"}) ./ (T_stats{idx, "ADC_avg"} - T_stats{GnoIPC, "ADC_avg"});
T_stats_normalized.DEL_std = 100*(T_stats.DEL_std - 0*T_stats{GnoIPC, "DEL_std"}) ./ (T_stats{idx, "DEL_avg"} - T_stats{GnoIPC, "DEL_avg"});



for i = 1:3
    G = Gs(:, i);
    %     x = xs(:, i);
    %     y = ys(:, i);
    patchcolor = patchcolors(i, :);
    
    xMean = T_stats{G, "ADC_avg"};
    xMeanNoIPC = T_stats{GnoIPC, "ADC_avg"};
    xMeanFullIPC = T_stats{idx, "ADC_avg"};
    xMean_rel = (xMean - xMeanNoIPC) ./ (xMeanFullIPC-xMeanNoIPC);
    
    yMean = T_stats{G, "DEL_avg"};
    yMeanNoIPC = T_stats{GnoIPC, "DEL_avg"};
    yMeanFullIPC = T_stats{idx, "DEL_avg"};
    yMean_rel = (yMean - yMeanNoIPC) ./ (yMeanFullIPC-yMeanNoIPC);
    
    [xMean_rel, isort] = sort(xMean_rel);
    yMean_rel = yMean_rel(isort);
    plot(ax, xMean_rel*100, yMean_rel*100, '-x', 'markersize', markersize, 'color', patchcolor)
    
    
    %     patch(ax, (x-xMeanNoIPC)./(xMeanFullIPC-xMeanNoIPC)*100, (y  - yMeanNoIPC) ./ (yMeanFullIPC-yMeanNoIPC)*100, ...
    %         1, 'FaceColor', patchcolor, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    drawStandardDeviationPatch(ax, T_stats_normalized, G, patchcolor, alpha, transparencies(i));
end

scatter(ax, 100, 100, 75, 'filled', 'k', 'linewidth', 1)
errorbar(ax, 0, 0, ...
    (T_stats{GnoIPC, "DEL_std"}) ./ (yMeanFullIPC-yMeanNoIPC)*100, ...
    (T_stats{GnoIPC, "DEL_std"}) ./ (yMeanFullIPC-yMeanNoIPC)*100, ...
    (T_stats{GnoIPC, "ADC_std"}) ./ (xMeanFullIPC-xMeanNoIPC)*100, ...
    (T_stats{GnoIPC, "ADC_std"}) ./ (xMeanFullIPC-xMeanNoIPC)*100, ...
    'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'Color', 'k')
% errorbar(ax, 100, 100, ...
%     T_stats{idx, "DEL_std"} ./ (yMeanFullIPC-yMeanNoIPC)*100, ...
%     T_stats{idx, "DEL_std"} ./ (yMeanFullIPC-yMeanNoIPC)*100, ...
%     T_stats{idx, "ADC_std"}./ (xMeanFullIPC-xMeanNoIPC)*100, ...
%     T_stats{idx, "ADC_std"}./ (xMeanFullIPC-xMeanNoIPC)*100, ...
%     'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'Color', 'k')



% Labels and stuff.
xlabel('Actuator duty cycle increase (%)')
ylabel('Damage equivalent load decrease (%)')
legend('$\ell^\infty$-IPC', '', '$\ell^2$-IPC', '', 'full IPC', '', '$w_c = 0.2$ rad/s', 'no IPC', 'location', 'southeast', 'interpreter', 'latex')
xlim([-5, 115])
ylim([-20, 140])
axis equal

postplot(f, fullfile(figDir, 'DELvsADC_turbulent_normalized.pdf'), 'fontsize', 8, 'fontname', 'NimbusRomNo9L', 'linewidth', 1);

% postplot(f, fullfile(pptDir, 'DELvsADC_turbulent_relative.svg'), 'fontsize', 15, 'fontname', 'NimbusRomNo9L', 'linewidth', 1.5);


%% Just no and full IPC as opening slide
[f, ax] = preplot('paperFormat', 'ppt', 'linefrac', 0.65, 'aspectRatio', 1.618, 'colororder', mycolors);

scatter(ax, T_stats{GnoIPC, "ADC_avg"}, T_stats{GnoIPC, "DEL_avg"}, ...
    'SizeData', 40, 'MarkerEdgeColor', 'k', 'linewidth', 1.5)
scatter(ax, T_stats{GfullIPC, "ADC_avg"}, T_stats{GfullIPC, "DEL_avg"}, ...
    'filled', 'SizeData', 40, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'linewidth', 1.5)


xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, 'Damage equivalent load (MNm)')
ylim(ax, [11, 15])
xlim(ax, [0, 0.24])
xticks(0:0.02:0.24)

legend('no IPC', 'full IPC')

postplot(f, fullfile(pptDir, 'noVsFullIPC.svg'), 'fontsize', 15, 'fontname', 'NimbusRomNo9L', 'linewidth', 1.5);
