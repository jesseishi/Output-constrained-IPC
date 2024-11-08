function y = radps2rpm(u)
% RADPS2RPM
arguments
    u = 1  % Allows this function to be called without arguments, facilitating y = radps2rpm * u in the command window use in Simulink as gain somewhere.
end
y = 60 / (2*pi) * u;
end

