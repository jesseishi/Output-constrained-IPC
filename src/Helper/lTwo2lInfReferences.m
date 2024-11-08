function refLInf = lTwo2lInfReferences(refL2,originalLoad_ty,makePlot)
%LTWO2LINFREFERENCES
arguments
    refL2 timeseries  % l2 reference.
    originalLoad_ty double
    makePlot logical  % Output a figure to check the results.
end

% Tilt-yaw load is [y, x]
My = originalLoad_ty(1);
Mx = originalLoad_ty(2);

if Mx < 0 || My < 0
    error("Not yet configured for negative original loads.")
end

% This was not super straightforward to find out. But if you plot it you
% can see why it is true.
if Mx > My
    refLInf = sqrt(refL2.Data.^2 - min(refL2.Data./sqrt(2), My).^2);
elseif My > Mx
    refLInf = sqrt(refL2.Data.^2 - min(refL2.Data./sqrt(2), Mx).^2);
else
    refLInf = refL2.Data / sqrt(2);
end
refLInf = timeseries(refLInf, refL2.Time);

if makePlot
    figure
    hold on
    axis('equal')

    % Draw the original load and a line to the origin.
    plot(Mx, My, 'kx', 'MarkerSize', 10)
    plot([0, Mx], [0, My], 'k')
    set(gca,'ColorOrderIndex',1)

    % Helper lines
    yline(My, 'k--')
    xline(Mx, 'k--')
    plot([-1000, 1000], [-1000, 1000], 'k--')

    % Draw the circles.
    theta = linspace(0, 2*pi, 100);
    rs = unique(refL2.Data)';
    for i = 1:length(rs)
        r = rs(i);
        plot(r * cos(theta), r * sin(theta))
    end

    % Draw the squares.
    set(gca,'ColorOrderIndex',1)
    ls = unique(refLInf.Data);
    for i = 1:length(ls)
        l = ls(i);

        x = [l, -l, -l, l, l];
        y = [l, l, -l, -l, l];
        plot(x, y, '--')
    end
end

end

