%% Setup.
restoredefaultpath;
clearvars; close all; clc;

% Relative location of .fst file.
fast.FAST_SFuncDir     = 'C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\openfast\build\bin';
fast.FAST_InputFile    = 'FullDOF_U15_RotatingWind_PLExp0p07_T900.fst';
fast.FAST_directory    = 'Results\InputFiles\FullDOF_U15_RotatingWind_PLExp0p07_T900';

% Add the path for the S-function, openFAST Matlab toolbox, and ROSCO
% toolbox.
addpath(fast.FAST_SFuncDir);
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\matlab-toolbox'))
addpath(genpath('C:\Users\jhummel\OneDrive - Delft University of Technology\Documenten\Projects\openFAST\ROSCO\Matlab_Toolbox\Utilities'))
addpath('preplot-postplot/src')

% Add other necessary paths.
addpath(genpath('src'))

%% Set the cache folder for Simulink and Matlab.
% This keeps our workspace cleaner.
cacheFolder = 'src/TempCache';
if ~isfolder(cacheFolder)
    mkdir(cacheFolder)
end
set_param(0, 'CacheFolder', cacheFolder)

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
[Ki, dMtdThetatLIN, dMydThetayLIN, dMdTheta] = tuneIPCcontrollers(windSpeed_mps, DT, wc, optimal_azimuth_offset_rad, false, false, true, R, F, true);
fprintf('Gain for all IPC controllers: Ki = %.4g for a wind speed of %.1f m/s.\n', Ki, windSpeed_mps)

% TODO: tune Ki_Mphi
% + investigate the sign. I changed it from neg to pos because with a neg
% gain I would end up at +/- 180 deg rotation from where I wanted.
Ki_Mphi = 1e-4;

% Original load estimator
w_originalLoadEstimator = 1;


%% Notch design.
[gmin, damp] = tuneNotchFilter(DT);

%% Make reference loads.
originalLoad_ty = [1812, 611];  % Needed for l2 to linf reference load conversion (to make comparison easier).
loadExactlyInTheCorner = min(originalLoad_ty) * sqrt(2);  % *sqrt(2) because we will convert this l2 load to linf later.
% steps = [5000, 1.5*loadExactlyInTheCorner, loadExactlyInTheCorner, 0.5*loadExactlyInTheCorner, 0];  % For l2 vs linf.
% steps = [5000, 1347, 792, 792/2, 0];  % This lines up even nicer.
steps = [5000, 0, 0, 0, 5000];  % For rotating wind.
% steps = [5000, 1750, 1500, 1250, 1000, 750, 500, 250, 100, 0];  % For realistic wind case
% steps = [5000, 4000, 3000, 2000, 1875. 1750, 1625, 1500, 1375, 1250, 1125];  % Some additional data
% steps = linspace(1800, 0, 10);
% steps = [10000, 500];

% Load for no IPC goes up to 9000 kNm in turbulent conditions.
% Takes about 20 minutes per simulation. So overnight can do about 24.
% steps = [10000, 5000, 1950, 1800, 1650, 1500, 1350, 1200, 1050, 900, 750, 600, 450, 300, 150, 0];  % For realistic wind data.
% steps = linspace(1e4, 0, 21);
[ref_loads, ref_load_names, TMaxs] = generateReferenceLoads(steps, "reference_type", "step", 'show_plot', true, 'TMax', 900, 'step_timestep', 40, 'step_starttime', 300);

ref_load = ref_loads{1};
TMax = TMaxs{1};

%% Simulate constrained IPC
close all

models = dir('src/Models/cIPC_l*originalLoad.slx');
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
        TTname = sprintf('%s_r%s_%s', model, ref_load_name, replace(fast.FAST_InputFile, '.fst', '.csv'));
        fname = fullfile(outDir, TTname);
        if isfile(fname)
            fprintf('Already have %s, skipping.\n', TTname); pause(0.2)
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
            fprintf('Simulating %s with ref load %s.\n', model, ref_load_name); pause(0.2);
            %             continue
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

beep

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
    TTname = sprintf('%s_r%s_%s', modelname, ref_load_name, replace(fast.FAST_InputFile, '.fst', '.csv'));
    fname = fullfile(outDir, TTname);
    if isfile(fname)
        fprintf('Already have %s, skipping.\n', TTname); pause(0.2)
        continue
    end
    
    fprintf('#############################\n')
    fprintf('Simulating fullIPC with Ki %g.\n', Ki); pause(0.1);
    fprintf('#############################\n')
    %     continue
    simOut = sim(model, [0, TMax]);
    
    % Save data.
    
    saveastimetable(fname, simOut, OutList, fast);
end

beep
