function plotPareto(ax, T, xSignal, ySignal, opts)
%PLOTPARETO Summary of this function goes here
%   Detailed explanation goes here
arguments
    ax  % Axes
    T table
    xSignal string
    ySignal string
    opts.randomVariable string = [];
    opts.optimizationVariable string = [];
    opts.lineStyle string = '-'
end

x = T.(xSignal);
y = T.(ySignal);

[x, isort] = sort(x);
y = y(isort);

plot(ax, x, y, opts.lineStyle)
xlabel(xSignal)
ylabel(ySignal)

end

