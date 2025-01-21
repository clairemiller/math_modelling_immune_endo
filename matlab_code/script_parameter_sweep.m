%% Output directory setup 
% (needed for all sections)
dirname = "../../data/parameter_sweep_local";
if ~exist(dirname, 'dir')
    mkdir(dirname);
end

% Set up the parallel processing
pool = gcp('nocreate'); % If no pool, do not create new one
n_procs = 6;
if isempty(pool)
    parpool(n_procs);
end


%% Setup for parameter sweep
%===============================================================================
% Sweep parameters
N_samples = 101;
omega_vals_log = linspace(-6,  -4,  N_samples);
beta1_vals_log = linspace(-7,  -5,  N_samples);
rho0_vals = [0.01, 0.1, 0.2];
nSims = N_samples*N_samples*length(rho0_vals);

% Prepolute a cell array of the parameter combinations for the solver with NaNs
parmsTable = table(...
                NaN(nSims, 1), NaN(nSims, 1), NaN(nSims, 1), NaN(nSims, 1), ...
                'VariableNames',{'id', 'log10omega', 'log10beta1', 'rho0'});
i = 0;
for rho0 = rho0_vals
    for omega_log = omega_vals_log
        for beta1_log = beta1_vals_log
            i = i+1;
            parmsTable(i,:) = table(i,omega_log, beta1_log, rho0, ...
                'VariableNames', {'id', 'log10omega', 'log10beta1', 'rho0'});
        end
    end
end
assert(i==nSims);

% Initial conditions (for constant influx)
varNames = {'M0', 'M1',      'M2',  'K0', 'KA','E0','EF','EA'};
y0       = [9e5;   1e2;     4.5e4;   5e5;  2e4;   0;   0;   0];

% Timing parameters
t0 = 0;
dt = 0.5;
tConstMinEval = 2e3;
tEndConst = 5e3;
tEndCyclic = 14*28;

% DE solver options
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);


%% Run parameter sweep
%===============================================================================

% Run the system for each parameter combination and store result in array
fprintf("Running %i simulations...\n",nSims);
compTimeStart = datetime("now");

% Prefill output container for parfor
output_table_constant = cell(nSims, 1);
output_table_cyclic = cell(nSims,1);

% Now loop over and store the simulation data (we do the analysis later)
parfor i = 1:nSims
    % Extract the correct parameter row
    parmMapi = struct( ...
        'rho0',parmsTable(i,:).rho0, ...
        'beta1',10^parmsTable(i,:).log10beta1, ...
        'omega',10^parmsTable(i,:).log10omega)

    % Solve constant ODE system at fixed intervals to calculate the median
    [tConst,yConst] = ode15s( ...
                @(t,y)fn_immune_ode_system(t, y, parmMapi, true), ...
                (t0:dt:tEndConst), y0', opts);
    median_yConst = median(yConst(tConst >= tConstMinEval,:));

    % Store the results of interest from the constant influx system
    metrics_constant = array2table(median_yConst,"VariableNames",varNames);
    [LowEA,~,~] = fn_determine_condition_satisfied(metrics_constant, parmMapi);
    metrics_constant.LowEA = LowEA;
    metrics_constant.id = parmsTable(i,:).id;
    output_table_constant{i} = metrics_constant;

    % Run the cyclic system at fixed intervals so we can calculate metrics
    new_y0 = median_yConst;
    [tCyclic,yCyclic] = ode15s( ...
                    @(t,y)fn_immune_ode_system(t, y, parmMapi, false), ...
                    (t0:dt:tEndCyclic), new_y0, opts);

    % Filter to time region of interest
    Y_filtered = yCyclic(tCyclic >= 10*28 & tCyclic < 13*28,:);
    
    % Compute metrics for each variable
    meanVals = mean(Y_filtered, 1);
    medianVals = median(Y_filtered, 1);
    finalVals = Y_filtered(end, :);
    
    % Create a table with metrics and add parameters
    metrics_table = array2table([meanVals; medianVals; finalVals], ...
                                "VariableNames",varNames);
    metrics_table.metric = ["mean"; "median"; "final"];
    metrics_table.id = repmat(parmsTable(i,:).id,3,1);
    output_table_cyclic{i} = metrics_table;

    % Somewhat see how we're tracking
    if (mod(i,500)==0)
        fprintf("i = %i\n", i)
    end
end
comptime = datetime("now")-compTimeStart


%% Save output data
%===============================================================================
disp("Writing data to file...");, 

% Concatenate the tables
table_constant = vertcat(output_table_constant{:});
table_cyclic = vertcat(output_table_cyclic{:});

% Write tables
writetable(parmsTable,fullfile(dirname, "parameters.csv"));
writetable(table_constant,fullfile(dirname, "constant_influx_summarydata.csv"));
writetable(table_cyclic,fullfile(dirname, "cyclic_influx_summarydata.csv"));

% Write default parameters to yaml
[muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
     KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
        rhoF, muE,omega,beta1,rho0] = default_parameters;
defaultParmsTable = table(muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0, ...
    deltaE,sigma,MC,KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK, ...
    gamma, rhoF, muE,omega,beta1,rho0);
fileID = fopen(fullfile(dirname,'default_parameters.yml'),'w');
for p = defaultParmsTable.Properties.VariableNames
    fprintf( fileID, "%s: %.3e\n", p{1}, defaultParmsTable.(p{1}) );
end
fclose(fileID);

% Write the run details as well
fileID2 = fopen(fullfile(dirname,'run_stats.yml'),'w');
fprintf( fileID2, "run_date_time: %s\n", compTimeStart);
[status, gitHash] = system('git rev-parse HEAD');
fprintf( fileID2, "git_hash: %s", gitHash);
[~, name] = system('uname -n');
fprintf( fileID2, "system_name: %s", name);
fprintf( fileID2, "n_procs: %i\n", n_procs);
fprintf( fileID2, "n_simulations: %i\n", nSims);
fprintf( fileID2, "comp_time: %s\n", comptime);
fclose(fileID2);
