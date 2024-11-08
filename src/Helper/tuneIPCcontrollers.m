function [Ki, dMtdThetat, dMydThetay, dMdTheta] = tuneIPCcontrollers( ...
    windSpeed_mps, DT, wc, azimuth_offset_rad, show_plot, ...
    useCPC, R, F, plotCPC)
% Tune IPC controllers based on linFiles.
% This controller targets DC signals and there is noise at 3kP. So at 3P
% and higher we should have good attenuation. 1P is about 10 rpm = 1 rad/s,
% so let's put the open-loop gain crossover a bit below that.
% On the other hand, I also use the notch filter to filter out the 3P
% signals so maybe a higher bandwidth is possible.
% Anyway, much faster than this with just an integrator is also not
% possible, but if I want to go faster I can omit wc and specify a phase
% margin.
% TODO: Take collective pitch control into account.
arguments
    windSpeed_mps double
    DT double
    wc double
    azimuth_offset_rad double = 0
    show_plot logical = false
    useCPC logical = false
    R = struct
    F = struct
    plotCPC logical = true
end

windSpeed_str = replace(sprintf('%.1f', windSpeed_mps), '.', 'p');
linFileNames = dir(sprintf('src/Lib/Linearizations/LIN_ws_%s_mps.*.lin', windSpeed_str));
linFileNames = fullfile({linFileNames.folder}, {linFileNames.name});
if isempty(linFileNames)
    error("Couldn't find linearizations for a wind speed of %.1f m/s", windSpeed_mps);
end
[MBC, matData, FAST_linData] = fx_mbc3_offset(linFileNames, azimuth_offset_rad);

% Take the average matrices to get an LTI system and remove the azimuth state (because
% it has no meaning in the nonrotating frame).
A = MBC.AvgA(2:end, 2:end);
B = MBC.AvgB(2:end, :);
C = MBC.AvgC(: , 2:end);
D = MBC.AvgD;
sys = ss(A, B, C, D, 'StateName', matData.DescStates(2:end), 'InputName', matData.DescCntrlInpt, 'OutputName', matData.DescOutput);

% Remove those AD User prop inputs.
iuKeep = ~contains(matData.DescCntrlInpt, 'AD User');
sys = sys(:, iuKeep);

% Add collective pitch control to this model before designing the
% individual pitch controller.
if useCPC
    iTheta_cpc = find(contains(matData.DescCntrlInpt, 'collective blade-pitch command'));
    iGenspeed = find(contains(matData.DescOutput, 'GenSpeed'));
    
    % Generator speed filter from ROSCO.
    GenSpeedFilter_d = tf(F.HSS.b, F.HSS.a, DT);
    GenSpeedFilter = d2c(GenSpeedFilter_d, 'tustin');
    GenSpeedFilter = GenSpeedFilter * rpm2radps;
    
    % Extract the operating condition.
    u_op = zeros(matData.NumInputs, matData.NAzimStep);
    for i = 1:matData.NAzimStep
        u_op(:, i) = cell2mat(FAST_linData(i).u_op);
    end
    Avgu_op = mean(u_op, 2);
    beta = Avgu_op(iTheta_cpc);
    if isfield(matData, 'Avguop')  % I've added this calculation to my local calculation of fx_mbc3_offset before I realized that it is also present in FAST_linData.
        beta2 = matData.Avguop(iTheta_cpc);
        assert(beta == beta2, sprintf("Both methods should yield the same operating condition. Got %.1f and %.1f rad.", beta, beta2));
    end
    
    % Gain scheduling from ROSCO
    PC_beta = R.PC_GS_angles;
    Kp_pc = R.PC_GS_KP;
    Ki_pc = R.PC_GS_KI;
    Kp      = interp1(PC_beta,Kp_pc, beta);
    Ki      = interp1(PC_beta,Ki_pc, beta);
    
    s = tf('s');
    C_CPC = Kp + Ki / s;
    
    plantWithController = sys;
    plantWithController(iGenspeed, iTheta_cpc) = C_CPC * plantWithController(iGenspeed, iTheta_cpc);
    sys_CPC_CL = feedback(plantWithController, GenSpeedFilter, iTheta_cpc, iGenspeed, -1);
    
    if plotCPC
        step(sys(iGenspeed, 1), sys_CPC_CL(iGenspeed, 1))  % Step response to a step in wind.
        title('Step response to a step in wind speed')
        ylabel('Generator speed (rad/s)')
    end
    
    % Put this in sys so that we use it to tune the IPC controllers.
    sys = sys_CPC_CL;
end


% Select the tilt and yaw moments and pitch as input and output.
% TODO: Select the indices programatically.
sys = sys(15:16, 5:6);

% Select one of the diagonal systems (doesn't matter which one) and design a controller
% for it in discrete time.
sys22 = sys(2, 2);
sys22m = minreal(sys22);
sys22md = c2d(sys22m, DT, 'zoh');  % Closest match to the openFAST model because it sees piecewise constant inputs from our controller, right?
[C_i, info] = pidtune(sys22md, 'I', wc);
if ~info.Stable
    error('Closed-loop is unstable, lower crossover frequency')
end
Ki = get(C_i, 'Ki');  % Used by Simulink.

% Also get the steady-state gain, which is used by the original load
% estimator.
dMtdThetat = dcgain(sys(1,1));
dMydThetay = dcgain(sys(2,2));
dMdTheta = mean([dMtdThetat, dMydThetay]);

if show_plot
    figure
    bodemag(sys)
    
    figure
    bode(sys22, sys22m, sys22md, sys22md*C_i)
    
    figure
    nyquist(sys22md*C_i)
end

end
