function [gmin, damp] = tuneNotchFilter(DT, show_plot)
arguments
    DT
    show_plot logical = false
end

% From trial and error:
gmin = db2mag(-40);
damp = 0.1;

% Frequency between 7.13355 - 7.13774 rpm = 0,74702 - 0,74746 rad/s at 10
% m/s wind speed without turbulence.
freqs = [0.74702, 0.74746];
freq = 3*mean(freqs);

N = tf([1, 2*gmin*damp*freq, freq^2], [1, 2*damp*freq, freq^2]);
Nd = c2d(N, DT, 'tustin');

if show_plot
    plotoptions = bodeoptions('cstprefs');
    plotoptions.MagUnits = 'abs';
    plotoptions.Grid = 'on';
    bode(N, Nd, plotoptions)
    legend
    xline(freq/3)
    xline(freq)
    xline(freqs)
end

end
