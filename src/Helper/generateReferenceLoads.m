function [y, y_name, TMaxs] = generateReferenceLoads(loads, opts)
% GENERATEREFERENCELOAD
% Generate one or multiple reference loads. It has two main usecases:
% 1) Generating multiple reference load timeseries all with its own load.
% 2) Generating a step reference that goes through all loads.
% 3) Both.
arguments
    loads (1,:) double
    opts.DT double = 0.005
    opts.TMax double = 300  % Only used for constant references.
    opts.reference_type string {mustBeMember(opts.reference_type, ["step", "constant", "both"])} = "step"
    opts.step_timestep double = 25
    opts.step_starttime double = 100
    opts.show_plot logical = false
end

switch opts.reference_type
    case "step"
        [y, y_name, TMaxs] = make_step_reference(loads, opts);
    case "constant"
        [y, y_name, TMaxs] = make_constant_reference(loads, opts);
    case "both"
        [y, y_name, TMaxs] = make_constant_reference(loads, opts);
        % Append step references to this.
        [y(end+1), y_name(end+1), TMaxs(end+1)] = make_step_reference(loads, opts);
    otherwise
        error("Reference type (%s) not recognised", opts.reference_type);
end

if opts.show_plot
    figure; hold on; grid on;
    for i = 1:length(y)
        plot(y{i})
    end
    legend(y_name)
end

end

function [y, y_name, TMaxs] = make_step_reference(loads, opts)
% Make the time array such that we can fit all the loads. NOTE: a user
% input on TMax will thus be ignored.
TLastStep = opts.step_starttime + opts.step_timestep * length(loads);
assert(opts.TMax >= TLastStep, "TMax truncates some of the steps, the last step ends at %.0f sec.", TLastStep)
t = 0:opts.DT:opts.TMax;

y = loads(1) * ones(size(t));
for i = 1:length(loads)
    y(t > opts.step_starttime + opts.step_timestep * (i-1)) = loads(i);
end
y = {timeseries(y, t, 'name', 'reference load step')};
y_name = {'step'};
TMaxs = {opts.TMax};
end

function [y, y_name, TMaxs] = make_constant_reference(loads, x)
% Use the default time length.
TMax = x.TMax;
t = 0:x.DT:TMax;

y = cell(size(loads));
y_name = cell(size(loads));
TMaxs = cell(size(loads));
for i = 1:length(loads)
    temp = loads(i) * ones(size(t));
    name = num2str(loads(i));
    y{i} = timeseries(temp, t, 'name', sprintf('reference load %s', name));
    y_name{i} = name;
    TMaxs{i} = TMax;
end
end