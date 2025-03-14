%% Setup.
restoredefaultpath;
clearvars; close all; clc;

% Set the cache folder for Simulink and Matlab.
% This keeps our workspace cleaner.
cacheFolder = 'src/TempCache';
if ~isfolder(cacheFolder)
    mkdir(cacheFolder)
end
set_param(0, 'CacheFolder', cacheFolder)

% Add the path for the S-function, openFAST Matlab toolbox, and ROSCO
% toolbox.
addpath(genpath('..\matlab-toolbox'))
addpath(genpath('..\ROSCO\Matlab_Toolbox\Utilities'))

% Add other necessary paths.
addpath(genpath('src'))


%% Loop through all input files.
inputFiles = [dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed100_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed101_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed102_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed103_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed104_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed105_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed106_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed107_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed108_T2100/*.fst");
    dir("Results/InputFiles/FullDOF_U15_TI8_PLExp0p07_seed109_T2100/*.fst")];

for iInputFile = 1:length(inputFiles)
    
    %% Define fast variable with input files.
    fast.FAST_SFuncDir     = '..\openfast-v3.5.0\build\bin';
    fast.FAST_InputFile    = inputFiles(iInputFile).name;
    fast.FAST_directory    = inputFiles(iInputFile).folder;
    addpath(fast.FAST_SFuncDir);
    fprintf('Selecting input file %s.\n', fast.FAST_InputFile)
    
    
    %% Additional setup.
    % Read FAST parameters.
    [FastParam, Cx] = ReadWrite_FAST(fast);
    
    % Get DT from the main input file.
    DT = FastParam.FP.Val{matches(FastParam.FP.Label,'DT')};
    
    % Get the generator efficiency and generator gain for region 2 control from
    % the ElastoDyn input file.
    GenEff = 1/100 * FastParam.SvDP.Val{matches(FastParam.SvDP.Label, 'GenEff')};
    VS_Rgn2K_Nmprad2ps2 = FastParam.SvDP.Val{matches(FastParam.SvDP.Label, 'VS_Rgn2K')} * radps2rpm^2;
    
    % Indexing.
    FAST_indexing = GetOutDataIndices();
    
    %% Load the collective pitch controller from ROSCO
    simu.dt = DT;
    [R, F] = load_ROSCO_params(FastParam, simu);
    
    
    %% Get the optimal azimuth offset.
    % windSpeed_mps = FastParam.IWP.Val{matches(FastParam.IWP.Label, 'HWindSpeed')};  % Doesn't work when using a .bts wind file.
    windSpeed_mps = 15.0;  % hard-coded!
    [optimal_azimuth_offset_rad, dMtdThetat, dMydThetay] = getOptimalAzimuthOffset(windSpeed_mps);
    
    % Convert to magnitude and flip the sign. The sign needs to be flipped
    % because it has a 180 deg phase offset.
    dMtdThetat = -db2mag(dMtdThetat);
    dMydThetay = -db2mag(dMydThetay);
    
    %% Tune the IPC controllers.
    wc = 0.2;  % gain crossover frequency, rad/s.
    [Ki, ~, ~, dMdTheta] = tuneIPCcontrollers(windSpeed_mps, DT, wc, optimal_azimuth_offset_rad, false, false);
    fprintf('Gain for all IPC controllers: Ki = %.4g for a wind speed of %.1f m/s.\n', Ki, windSpeed_mps)
    
    % TODO: tune Ki_Mphi
    % + investigate the sign. I changed it from neg to pos because with a neg
    % gain I would end up at +/- 180 deg rotation from where I wanted.
    Ki_Mphi = 1e-4;
    
    % Original load estimator
    w_originalLoadEstimator = 2;
    
    
    %% Notch design.
    [gmin, damp] = tuneNotchFilter(DT);
    
    %% Make reference loads.
    % U10 original load (545, 1435)
    % U15 original load (1350, 500) (tilt, yaw) with peaks to (1800, 700)
    originalLoad_ty = [1812, 611];  % Needed for l2 to linf reference load conversion (to make comparison easier).
    % steps = [5000, 1000, 500, 0];  % For l2 vs linf.
    % steps = [5000, 0, 0, 0, 5000];  % For rotating wind.
    % steps = [5000, 1750, 1500, 1250, 1000, 750, 500, 250, 100, 0];  % For realistic wind case
    % steps = [5000, 4000, 3000, 2000, 1875. 1750, 1625, 1500, 1375, 1250, 1125];  % Some additional data
    
    % Load for no IPC goes up to 9000 kNm in turbulent conditions.
    % Takes about 20 minutes per simulation. So overnight can do about 24.
    steps = [10000, 4000:-500:0];
    %     steps = [10000, 3000, 2500, 2000, 1500, 0];
    [ref_loads, ref_load_names, TMaxs] = generateReferenceLoads(steps, "reference_type", "constant", 'show_plot', true, 'TMax', 2100, 'step_timestep', 40, 'step_starttime', 100);
    
    ref_load = ref_loads{1};
    TMax = TMaxs{1};
    
    %% Simulate constrained IPC
    close all
    
    models = dir('src/Models/cIPC_leakyIntegrator.slx');
    errors = [];
    
    for i = 1:length(models)
        for j = 1:length(ref_loads)
            model = replace(models(i).name, '.slx', '');
            ref_load = ref_loads{j};
            ref_load_name = ref_load_names{j};
            TMax = TMaxs{j};
            
            % Skip the run if we've already simulated it.
            % Save the logged signals.
            temp = replace(fast.FAST_directory, 'InputFiles', 'data');
            outDir = sprintf('%s/ref%s', temp, ref_load_name);
            if ~isfolder(outDir)
                mkdir(outDir)
            end
            TTname = sprintf('%s_r%s_w%.1f_wo%.1f_%s', model, ref_load_name, wc, w_originalLoadEstimator, replace(fast.FAST_InputFile, '.fst', '.csv'));
            fname = fullfile(outDir, TTname);
            if isfile(fname)
                fprintf('Already have %s, skipping.\n', TTname); pause(0.02)
                continue
            end
            
            
            % If we use an l-inf based controller, we need to adjust the
            % reference, so that the l-2 norm that it produces remains the
            % same. NOTE: only works for wind files where the shear and veer
            % don't change, such as the rotating wind case, and strictly
            % speaking the turbulent case.
            convertl2tolinf = true;
            if convertl2tolinf
                if contains(model, 'inf')
                    ref_load = lTwo2lInfReferences(ref_load, originalLoad_ty, false);
                end
            else
                warning('Not converting l2 to linf reference.'); pause(1);
            end
            
            
            % Try simulating, if there's an error save it but continue to the
            % next run.
            try
                fprintf('Simulating %s with ref load %s.\n', model, ref_load_name); pause(0.02);
                %                 continue
                simOut = sim(model, [0, TMax]);
                
                
                % Save data.
                saveastimetable(fname, simOut, OutList, fast);
            catch e
                fprintf('Error while simulating %s with ref load %s.\n', model, ref_load_name);
                disp(e)
                fprintf('But continuing now with the next one.\n')
                beep
                errors = [errors, e];
            end
        end
    end
    
    %     beep
    
    %% Also do a full and no IPC run for this input file.
    % For noIPC we also use the fullIPC model, but just set the gains to 0.
    close all;
    
    ref_load_name = 'none';
    Kis = [0, Ki];
    
    for i = 1:length(Kis)
        Ki = Kis(i);
        
        model = 'IPC_fullIPC';
        if Ki == 0
            modelname = 'IPC_noIPC';
        else
            modelname = 'IPC_fullIPC';
        end
        
        TMax = max(cell2mat(TMaxs));
        
        % Skip the run if we've already simulated it.
        % Save the logged signals.
        temp = replace(fast.FAST_directory, 'InputFiles', 'data');
        outDir = sprintf('%s/ref%s', temp, ref_load_name);
        if ~isfolder(outDir)
            mkdir(outDir)
        end
        TTname = sprintf('%s_r%s_w%.1f_wo%.1f_%s', modelname, ref_load_name, wc, w_originalLoadEstimator, replace(fast.FAST_InputFile, '.fst', '.csv'));
        fname = fullfile(outDir, TTname);
        if isfile(fname)
            fprintf('Already have %s, skipping.\n', TTname); pause(0.02)
            continue
        end
        
        fprintf('#############################\n')
        fprintf('Simulating fullIPC with Ki %g.\n', Ki); pause(0.01);
        fprintf('#############################\n')
        %         continue
        simOut = sim(model, [0, TMax]);
        
        % Save data.
        
        saveastimetable(fname, simOut, OutList, fast);
    end
    
    %     beep
end

beep