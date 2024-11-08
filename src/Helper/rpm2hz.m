function y = rpm2hz(u)
arguments
    u = 1  % Allows this function to be called without arguments, facilitating y = rpm2hz * u in the command window use in Simulink as gain somewhere.
end
y = u / 60;
end
