function y = rpm2radps(u)
%RPM2RADPS
arguments
    u = 1  % Allows this function to be called without arguments, facilitating y = rpm2radps * u in the command window use in Simulink as gain somewhere.
end
y = 2*pi/60 * u;
end

