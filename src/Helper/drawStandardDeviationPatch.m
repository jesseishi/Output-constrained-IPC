function drawStandardDeviationPatch(ax, T_stats, G, basecolor, alpha, brightness)
xLinf = T_stats{G, "ADC_avg"} - T_stats{G, "ADC_std"};
xLinf = [xLinf; T_stats{G, "ADC_avg"}(end) + T_stats{G, "ADC_std"}(end)];
xLinf = [xLinf; flipud(T_stats{G, "ADC_avg"}) + flipud(T_stats{G, "ADC_std"})];
xLinf = [xLinf; xLinf(1); xLinf(1)];

yLinf = T_stats{G, "DEL_avg"} - T_stats{G, "DEL_std"};
yLinf = [yLinf; T_stats{G, "DEL_avg"}(end) - T_stats{G, "DEL_std"}(end)];
yLinf = [yLinf; flipud(T_stats{G, "DEL_avg"}) + flipud(T_stats{G, "DEL_std"})];
yLinf = [yLinf; T_stats{G, "DEL_avg"}(1) + T_stats{G, "DEL_std"}(1)];
yLinf = [yLinf; yLinf(1)];

basecolor = rgb2hsv(basecolor);
basecolor(3) = brightness;
basecolor = hsv2rgb(basecolor);
patch(ax, xLinf, yLinf, 1, 'FaceColor', basecolor, 'EdgeColor', 'none', 'FaceAlpha', alpha);
end