%% Output directory setup 
dirname = "../output/parameter_sweep_rhoF";
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
rhoF_vals = [0.01, 0.1, 0.2, 0.5];
nSims = N_samples*N_samples*length(rhoF_vals);

% Prepolute a cell array of the parameter combinations for the solver with NaNs
parmsTable = table(...
                NaN(nSims, 1), NaN(nSims, 1), NaN(nSims, 1), NaN(nSims, 1), ...
                'VariableNames',{'id', 'log10omega', 'log10beta1', 'rhoF'});
i = 0;
for rhoF = rhoF_vals
    for omega_log = omega_vals_log
        for beta1_log = beta1_vals_log
            i = i+1;
            parmsTable(i,:) = table(i,omega_log, beta1_log, rhoF, ...
                'VariableNames', {'id', 'log10omega', 'log10beta1', 'rhoF'});
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
        'rhoF',parmsTable(i,:).rhoF, ...
        'beta1',10^parmsTable(i,:).log10beta1, ...
        'omega',10^parmsTable(i,:).log10omega)

    % Solve constant ODE system at fixed intervals to calculate the median
    [tConst,yConst] = ode15s( ...
                @(t,y)immune_ode_system_rhoF(t, y, parmMapi, true), ...
                (t0:dt:tEndConst), y0', opts);
    median_yConst = median(yConst(tConst >= tConstMinEval,:));

    % Store the results of interest from the constant influx system
    metrics_constant = array2table(median_yConst,"VariableNames",varNames);
    [LowEA,~,~] = determine_condition_satisfied_rhoF(metrics_constant, parmMapi);
    metrics_constant.LowEA = LowEA;
    metrics_constant.id = parmsTable(i,:).id;
    output_table_constant{i} = metrics_constant;

    % Run the cyclic system at fixed intervals so we can calculate metrics
    new_y0 = median_yConst;
    [tCyclic,yCyclic] = ode15s( ...
                    @(t,y)immune_ode_system_rhoF(t, y, parmMapi, false), ...
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

%% Custom function for immune model to input rhoF
% t: current time, Y: y, p: params (type struct)
function dYdt = immune_ode_system_rhoF(t, Y, p, constMuE)
% Y order: [M0, M1, M2, K0, KA, E0, EF, EA]

    % Parameters --------------------------------------------
    [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
         KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma,rhoF, ...
            muE,omega,beta1,rho0] = default_parameters;
    
    % Input parameters
    omega = p.omega;
    beta1 = p.beta1;
    rhoF = p.rhoF; 

    % Calculate muE
    if (~constMuE)
        muE = (2*10^4).*(sin(pi*(t+12)/28)).^160;
    end

    % Calculate the system of DEs ----------------------------
    % Macrophage
    dM0dt = muM+etaM*(1-(Y(1)+Y(2)+Y(3))/MC)*(Y(1)+Y(2)+Y(3))-beta1*(Y(7)+Y(8))*Y(1)-thetaM*(Y(5)/(CKA+Y(5)))*Y(1)-beta2*Y(1)-deltaM*Y(1);
    dM1dt = beta1*(Y(7)+Y(8))*Y(1)+thetaM*(Y(5)/(CKA+Y(5)))*Y(1)+beta21*Y(3)-beta12*Y(2)-deltaM*Y(2);
    dM2dt = beta2*Y(1)+beta12*Y(2)-beta21*Y(3)-deltaM*Y(3);
    % Natural killer
    dK0dt = muK-thetaK*Y(4)*(Y(2)/(CM1+Y(2)))-deltaK*Y(4);
    dKAdt = thetaK*Y(4)*(Y(2)/(CM1+Y(2)))+etaK*(1-(Y(5)+Y(4))/KC)*Y(5)-sigma*(Y(7)+Y(8))*Y(5)-deltaK*Y(5);
    % Endometrial
    dE0dt = muE-deltaE0*(rho0*Y(6)+(1-rho0)*Y(6));
    dEFdt = deltaE0*rho0*Y(6)-omega*(gamma*Y(5)+(1-gamma)*Y(2))*Y(7)-rhoF*Y(7)-deltaE*Y(7);
    dEAdt = rhoF*Y(7)+etaE*(1-Y(8)/EC)*(Y(3)/(CM2+Y(3)))*Y(8)-omega*(gamma*Y(5)+(1-gamma)*Y(2))*Y(8)-deltaE*Y(8);

    % Concatenate into column array for return
    dYdt = [dM0dt; dM1dt; dM2dt; ...
            dK0dt; dKAdt; ...
            dE0dt; dEFdt; dEAdt];
end

%% Function for determining the EA condition
function [lowEASystem, growth_rate, clear_rate] = ...
                                    determine_condition_satisfied_rhoF(T,parms)
        
    % Read in the default parameters
    [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
              KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
                 rhoF, muE,omega,beta1,rho0] = default_parameters;
    
    % Overwrite the omega, beta1, and rho0 values
    omega = parms.omega;
    beta1 = parms.beta1;
    rhoF = parms.rhoF;
    
    % Calculate the two sides of the equation
    M2_hat = T.M2./(CM2+T.M2);
    growth_rate = etaE*M2_hat;
    clear_rate = omega * (gamma*T.KA + (1-gamma)*T.M1) + deltaE;
    
    % Determine if the condition is true
    lowEASystem = (growth_rate < clear_rate);
end
