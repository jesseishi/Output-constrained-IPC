%% Different norms
mycolors = [
    38,  87, 109;
         53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;
[f, ax] = preplot('paperFormat', 'ppt', 'lineFrac', 0.5, 'aspectRatio', 1, ...
    'interpreter', 'latex', 'colororder', mycolors);

% L1 norm
x = [1, 0, -1, 0, 1];
y = [0, 1, 0, -1, 0];
plot(x, y);

% L2 norm
phi = linspace(0, 2*pi, 1e3);
r = 1;
x = r*cos(phi);
y = r*sin(phi);
plot(x, y);

% Linf norm
x = [1, -1, -1, 1, 1];
y = [1, 1, -1, -1, 1];
plot(x, y);

% Labels and stuff
xlim([-1.5, 1.5])
axis(ax, 'equal')
xlabel('x (-)')
ylabel('y (-)')
legend('$\ell^1$-norm', '$\ell^2$-norm', '$\ell^\infty$-norm')

postplot(f, 'Results/figures/NAWEA/norms.svg', 'fontsize', 15, 'linewidth', 1.5)



%% Plot the transition from no IPC to full IPC.
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% nameFilter = 'Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refstep/*.csv';
[TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(...
    "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnoToFullIpc/cIPC_l*originalLoad*.csv", ...
    "TStart", seconds(300), 'TEnd', seconds(380));

figDir = 'Results/figures/NAWEA';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

%% Generic check.
plotTTsOverview(TTs)

%% Plot in the rotating frame
% There's only 1 timetable
TT = TTs{1};
[f, axs] = preplot(3, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors, 'Padding', 'compact');

plot(axs(1), TT(S, :), "FlappingMoment1_MNm")
plot(axs(2), TT(S, :), "FlappingMoment1_1Pfilt_MNm")
% plot(axs(2), TT(S, :).Time, TT(S, :).ref_load./1e3, 'k--')
plot(axs(3), TT(S, :), "BldPitch1_deg")

% Set limits and ticks.
ylim(axs(2), [-2, 2])

% Set the labels.
ylabel(axs(1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2), {'1P flapping', 'moment (MNm)'})
ylabel(axs(3), 'Pitch angle (deg)')
xlabel(axs(3), 'Time (s)')
title(axs(1), 'Rotating frame')

% legend(axs(2), '', 'reference')

postplot(f, fullfile(figDir, 'timeResponseRotatingFrame.svg'), ...
    'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5);


%% Plot in the nonrotating frame
[f, axs] = preplot(4, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors, 'Padding', 'compact');

plot(axs(1), TT(S, :), "M_tilt_filt_MNm")
% plot(axs(1), TT(S, :), "ref_load_tilt", 'color', 'k', 'linestyle', '--')
plot(axs(2), TT(S, :), "M_yaw_filt_MNm")
% plot(axs(2), TT(S, :), "ref_load_yaw", 'color', 'k', 'linestyle', '--')
plot(axs(3), TT(S, :), "theta_tilt_deg")
plot(axs(4), TT(S, :), "theta_yaw_deg")

% Set limits and ticks.
ylim(axs(1), [0, 2])

% Set the labels.
ylabel(axs(1), {'Tilt moment', '(MNm)'})
ylabel(axs(2), {'Yaw moment', '(MNm)'})
ylabel(axs(3), {'Tilt pitch', '(deg)'})
ylabel(axs(4), {'Yaw pitch', '(deg)'})
xlabel(axs(4), 'Time (s)')
title(axs(1), 'Nonrotating frame')

% legend(axs(1), '', 'reference')

% Share some y limits.
ylim(axs(1), [-0.010, 2])
linkprop([axs(1), axs(2)], 'YLim')
linkprop([axs(3), axs(4)], 'YLim')

postplot(f, fullfile(figDir, 'timeResponseNonrotatingFrame.svg'), ...
    'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5);



%% Laminar results
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% nameFilter = 'Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refstep/*.csv';
[TTs, Ts_fft, T, S, N, metadata] = loadandprocess2(...
    ["Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refstep/cIPC_l*originalLoad*.csv", ...
     "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/*.csv"], ...
    "TStart", seconds(300), 'TEnd', seconds(425), 'assertEqualLength', false);

figDir = 'Results/figures/NAWEA';
mycolors = [
    38,  87, 109;
    %      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;
model_labels = metadata(2, :);

iLinf = find(matches(model_labels, 'linforiginalLoad'));
iL2 = find(matches(model_labels, 'l2originalLoad'));
iOriginalLoad = [iLinf, iL2];

%% Generic check if the simulations look OK.
plotTTsOverview(TTs, 'legendLabels', model_labels)


%% Rotating frame for Linf
TT = TTs{iLinf};  % Select the Linf controller.
[f, axs] = preplot(3, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors, 'Padding', 'compact');

plot(axs(1), TT(S, :), "FlappingMoment1_MNm")
plot(axs(2), TT(S, :), "FlappingMoment1_1Pfilt_MNm")
plot(axs(2), TT(S, :).Time, TT(S, :).ref_load./1e3, 'k--')
plot(axs(3), TT(S, :), "BldPitch1_deg")

% Set limits and ticks.
ylim(axs(2), [-2, 2])

% Set the labels.
ylabel(axs(1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2), {'1P flapping', 'moment (MNm)'})
ylabel(axs(3), 'Pitch angle (deg)')
xlabel(axs(3), 'Time (s)')
title(axs(1), 'Rotating frame')

legend(axs(2), '', 'reference')

postplot(f, fullfile(figDir, 'timeResponseRotatingFrameLinf.svg'), ...
    'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5);


%% On the tilt-yaw plane for Linf
[f, axs] = preplot(2, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'colororder', mycolors, 'aspectRatio', 0.73);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

for i = 1:2
    set(axs(i), 'colororderindex', 1);
end

plot(axs(1), TT(S, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm")
plot(axs(2), TT(S, :), "theta_yaw_deg", "theta_tilt_deg")


ylabel(axs(1), 'Tilt moment (MNm)')
xlabel(axs(1), 'Yaw moment (MNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'Yaw pitch (deg)')
title(axs(1), 'Tilt-yaw plane')

% Let's also draw the reference loads in there.
linfrefs = unique(TTs{find(contains(model_labels, "linf"), 1)}.ref_load);
for i = 1:length(linfrefs)
    r = linfrefs(i);
    r = r./1e3;
    % Draw a box.
    x = [r, -r, -r, r, r];
    y = [r, r, -r, -r, r];
    plot(axs(1), x, y, 'k--')
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-1, 2])
ylim(axs(1), [-0.100, 2])
xlim(axs(2), [-0.2, 0.4])
ylim(axs(2), [0, 0.4]);

% legend(['no IPC'; 'full IPC'; model_labels(originalLoadIndices)], 'location', 'southeast')
% lg = legend(axs(1), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'northeast', 'interpreter', 'latex');
% lg.Layout.Tile = 'East';

% This is also not idea.
plot(axs(2), -10:-9, 0:1, '--k')  % Phantom line.
legend(axs(2), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'ty_response.svg'), ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)

%% Only the load
[f, ax] = preplot('lineFrac', 0.45, 'paperFormat', 'ppt', 'colororder', mycolors, 'aspectRatio', 1.35);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, TTs{inoIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(ax, TTs{ifullIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');


set(ax, 'colororderindex', 1);


plot(ax, TT(S, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm")


ylabel(ax, 'Tilt moment (MNm)')
xlabel(ax, 'Yaw moment (MNm)')
title(ax, 'Tilt-yaw plane')

% Let's also draw the reference loads in there.
linfrefs = unique(TTs{find(contains(model_labels, "linf"), 1)}.ref_load);
for i = 1:length(linfrefs)
    r = linfrefs(i);
    r = r./1e3;
    % Draw a box.
    x = [r, -r, -r, r, r];
    y = [r, r, -r, -r, r];
    plot(ax, x, y, 'k--')
end

axis(ax, 'equal')
xlim(ax, [-1, 2])
ylim(ax, [-0.100, 2])

% legend(['no IPC'; 'full IPC'; model_labels(originalLoadIndices)], 'location', 'southeast')
% lg = legend(axs(1), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'northeast', 'interpreter', 'latex');
% lg.Layout.Tile = 'East';

% This is also not ideal.
% plot(ax, -10:-9, 0:1, '--k')  % Phantom line.
legend(ax, 'no IPC', 'full IPC', '$\ell^\infty$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'ty_response_theta.svg'), ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)


%% Show on top of each other.
[f, axs] = preplot(3, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors, 'Padding', 'compact');

for i = 1:2
    j = iOriginalLoad(i);
    TT = TTs{j};
    plot(axs(1), TT(S, :), "FlappingMoment1_MNm")
    plot(axs(2), TT(S, :), "FlappingMoment1_1Pfilt_MNm")
    plot(axs(3), TT(S, :), "BldPitch1_deg")
end

plot(axs(2), TT(S, :).Time, TT(S, :).ref_load./1e3, 'k--')

% Set limits and ticks.
ylim(axs(2), [-2, 2])

% Set the labels.
ylabel(axs(1), {'Flapping', 'moment (MNm)'})
ylabel(axs(2), {'1P flapping', 'moment (MNm)'})
ylabel(axs(3), 'Pitch angle (deg)')
xlabel(axs(3), 'Time (s)')
title(axs(1), 'Rotating frame')

legend(axs(2), '', '', 'reference')

legend(axs(3), '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'timeResponseRotatingFrameBoth.svg'), ...
    'linkaxes', 'x', 'sharex', true, 'removeTimetableUnit', true, ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5);


[f, axs] = preplot(2, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'colororder', mycolors, 'aspectRatio', 0.73);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(axs(1), TTs{inoIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(axs(2), TTs{inoIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(axs(1), TTs{ifullIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
scatter(axs(2), TTs{ifullIPC}(end, :), "theta_yaw_deg", "theta_tilt_deg", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

for i = 1:2
    set(axs(i), 'colororderindex', 1);
end
for i = 1:2
    j = iOriginalLoad(i);
    TT = TTs{j};
    plot(axs(1), TT(S, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm")
    plot(axs(2), TT(S, :), "theta_yaw_deg", "theta_tilt_deg")
end

ylabel(axs(1), 'Tilt moment (MNm)')
xlabel(axs(1), 'Yaw moment (MNm)')
ylabel(axs(2), 'Tilt pitch (deg)')
xlabel(axs(2), 'Yaw pitch (deg)')
title(axs(1), 'Nonrotating frame')

% Let's also draw the reference loads in there.
l2refs = unique(TTs{find(contains(model_labels, "l2"), 1)}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i);
    r = r./1e3;
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(axs(1), r*cos(theta), r*sin(theta), 'k--')
end
linfrefs = unique(TTs{find(contains(model_labels, "linf"), 1)}.ref_load);
for i = 1:length(linfrefs)
    r = linfrefs(i);
    r = r./1e3;
    % Draw a box.
    x = [r, -r, -r, r, r];
    y = [r, r, -r, -r, r];
    plot(axs(1), x, y, 'k--')
end

axis(axs(1), 'equal')
axis(axs(2), 'equal')
xlim(axs(1), [-1, 2])
ylim(axs(1), [-0.100, 2])
xlim(axs(2), [-0.2, 0.4])
ylim(axs(2), [0, 0.4]);

% legend(['no IPC'; 'full IPC'; model_labels(originalLoadIndices)], 'location', 'southeast')
% lg = legend(axs(1), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'northeast', 'interpreter', 'latex');
% lg.Layout.Tile = 'East';

% This is also not idea.
plot(axs(2), -10:-9, 0:1, '--k')  % Phantom line.
legend(axs(2), 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'ty_response_both.svg'), ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)




% And again only theta
[f, ax] = preplot('lineFrac', 0.45, 'paperFormat', 'ppt', 'colororder', mycolors, 'aspectRatio', 1.35);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, TTs{inoIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(ax, TTs{ifullIPC}(end, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

set(ax, 'colororderindex', 1);

for i = 1:2
    j = iOriginalLoad(i);
    TT = TTs{j};
    plot(ax, TT(S, :), "M_yaw_filt_MNm", "M_tilt_filt_MNm")
end

ylabel(ax, 'Tilt moment (MNm)')
xlabel(ax, 'Yaw moment (MNm)')
title(ax, 'Tilt-yaw plane')

% Let's also draw the reference loads in there.
l2refs = unique(TTs{find(contains(model_labels, "l2"), 1)}.ref_load);
for i = 1:length(l2refs)
    r = l2refs(i);
    r = r./1e3;
    % Draw a circle;
    theta = linspace(0, 2*pi, 100);
    plot(ax, r*cos(theta), r*sin(theta), 'k--')
end
linfrefs = unique(TTs{find(contains(model_labels, "linf"), 1)}.ref_load);
for i = 1:length(linfrefs)
    r = linfrefs(i);
    r = r./1e3;
    % Draw a box.
    x = [r, -r, -r, r, r];
    y = [r, r, -r, -r, r];
    plot(ax, x, y, 'k--')
end

axis(ax, 'equal')
xlim(ax, [-1, 2])
ylim(ax, [-0.100, 2])

% This is also not idea.
legend(ax, 'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC', 'reference', 'location', 'southeast', 'interpreter', 'latex');


postplot(f, fullfile(figDir, 'ty_response_both_theta.svg'), ...
    'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)




%% Laminar trade-off
% Reset and load the appropriate things.
clearvars; close all; clc;

addpath(genpath('src'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath('preplot-postplot/src')

% We need all the constant references and no/full IPC.
nameFilter1 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/ref*0/cIPC_l*originalLoad*.csv";  % This feels a bit like a cheat. I should have a better way to index constant references.
nameFilter2 = "Results/data/FullDOF_U15_TI0_PLExp0p07_seed0_T900/refnone/*.csv";
[~, ~, T, ~, ~, metadata] = loadandprocess2([nameFilter1, nameFilter2], ...
    'assertEqualLength', false, 'TStart', seconds(300), 'TEnd', seconds(400), 'saveTTs', false, 'saveTsfft', false);
% Not all simulations were the same length -> Only use the part from 300 to
% 400 seconds. This is long enough since we're doing laminar flow.

figDir = 'Results/figures/NAWEA';
mycolors = [
     38,  87, 109;
%      53, 143, 161;
    109, 194, 155;
    188, 204, 151;
    235, 215, 155] ./ 255;

model_labels = metadata(2, :);

beep;

%% Generic check if the simulations look OK.
% plotTTsOverview(TTs)

%% Mean M vs mean theta
[f, ax] = preplot(1, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, T(inoIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(ax, T(ifullIPC, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

set(ax, 'ColorOrderIndex', 1)

unique_model_labels = unique(model_labels);
unique_cIPC_model_labels = unique_model_labels(~contains(unique_model_labels, 'IPC'));
unique_cIPC_model_labels = unique_cIPC_model_labels([2, 1]);  % Change the order so that linf comes first.
for i = 1:length(unique_cIPC_model_labels)
    model_label = unique_cIPC_model_labels(i);
    idx = matches(model_labels, model_label);

    plotPareto(ax, T(idx, :), "theta_ty_mag_deg_mean", "M_ty_mag_MNm_mean")
%     disp(idx)
end

xlabel(ax, 'Pitch magnitude (deg)')
ylabel(ax, 'Moment magnitude (MNm)')

legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'meanMvsT_laminar.svg'), 'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)

%% DEL vs ADC
[f, ax] = preplot(1, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, T(inoIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(ax, T(ifullIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

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
ylabel(ax, 'Damage equivalent load (MNm)')

legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_laminar.svg'), 'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)


%% Also calculate percentage-wise the gains.
iInf = matches(model_labels, 'linforiginalLoad');
iL2 = matches(model_labels, 'l2originalLoad');

xLInf = T{iInf, "theta_ty_mag_deg_mean"};
yLInf = T{iInf, "M_ty_mag_MNm_mean"};
xL2 = T{iL2, "theta_ty_mag_deg_mean"};
yL2 = T{iL2, "M_ty_mag_MNm_mean"};

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


%% Plot for the opening, just no and full IPC.
[f, ax] = preplot(1, 1, 'lineFrac', 0.45, 'paperFormat', 'ppt', 'aspectRatio', 1.1, 'colororder', mycolors);

% Add no and full IPC.
inoIPC = matches(model_labels, 'noIPC');
ifullIPC = matches(model_labels, 'fullIPC');
scatter(ax, T(inoIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'SizeData', 120, 'MarkerEdgeColor', 'k', 'lineWidth', 2);
scatter(ax, T(ifullIPC, :), "BldPitch1_deg_DC", "FlappingMoment1_MNm_DEL", 'filled', 'SizeData', 120, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');


xlabel(ax, 'Actuator duty cycle (-)')
ylabel(ax, 'Damage equivalent load (MNm)')

legend({'no IPC', 'full IPC', '$\ell^\infty$-IPC', '$\ell^2$-IPC'}, 'location', 'northeast', 'interpreter', 'latex')
postplot(f, fullfile(figDir, 'DELvsADC_laminar_conventional.svg'), 'fontname', 'NimbusRomNo9L', 'fontsize', 15, 'linewidth', 1.5)




%% Plot for turbulence
% See Plot_WES because that was 99% code overlap.















