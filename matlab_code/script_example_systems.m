%% Visualise an example
parms = struct("rho0", 0.01,  ...
               "omega", 10^-5, ...
               "beta1", 10^-6 )

% Solve
[yCyclic, yConst, y0Cyclic] = solve_system(parms);
plot_system_response( yConst )
plot_system_response( yCyclic )

plot(yConst.t, yConst.EA);
set(gca, 'YScale', 'log');
xlim([0 1000]);


% Print whether the condition is satisfied
T = median(yConst(yConst.t>1e3,:))
[lowEASystem, growth_rate, clear_rate] = fn_determine_condition_satisfied(T,parms)

% Print the scaled values for auto
[yConst(end,1:5)./1e6 yConst(end,6:8)./1e9]



%% Default parameters simulation
%===============================================================================
% Default parameters
defParms = struct();
defParms.rho0 = 0.1;
defParms.omega = 1e-5;
defParms.beta1 = 1e-6;

% Solve
[yCyclic, yConst, y0Cyclic] = solve_system(defParms);


% Plot system responses --------------------------------------------------------
close all
% Constant influx system
% Choose a time period
tmin = max(yConst.t) - 100;
iConst = (yConst.t >= tmin);
% Plot
plot_system_response( yConst )
sgtitle("Constant influx")

% Cyclic influx system
% Choose a time period
tmin = max(yCyclic.t) - 3*28;
iCyclic = (yCyclic.t >= tmin);
% Plot
plot_system_response( yCyclic )
sgtitle("Cyclic influx")


% Write (cyclic) result --------------------------------------------------------
folder = '../../data/default_system';
% Check if the folder exists and, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end
writetable(yCyclic, fullfile(folder,'timeseries.csv'));
write_sim_info(fullfile(folder,'model_info.yml'), ...
                    defParms, y0Cyclic, "defaultparms-cyclic");


%% No attachment example
parms = struct("rho0", 0.1,  ...
               "omega", 10^-4, ...
               "beta1", 10^-5 )

% Solve
[yCyclic, yConst, y0Cyclic] = solve_system(parms);
plot_system_response( yConst )
plot_system_response( yCyclic )


% Save
folder = '../../data/no_attachment_example';
% Check if the folder exists and, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end
writetable(yCyclic, fullfile(folder,'timeseries.csv'));
write_sim_info(fullfile(folder,'model_info.yml'), ...
                    parms, y0Cyclic, "defaultparms-cyclic");


%% Varying rho0 examples
%===============================================================================
% Parameters
rho0_vals = [0.0, 0.005, 0.01, 0.05, 0.1, 0.2];
parms = get_default_parms_struct();
parms = struct( 'omega',1e-5,  'beta1', 1e-6 );

% Output containers for the results
yOut = cell(length(rho0_vals),1);
y0_arr = zeros(length(rho0_vals),8);

% Loop over the rho0 values and store the timeseries for the cyclic
for i = 1:length(rho0_vals)
    parms.rho0 = rho0_vals(i);
    [yCyclic, ~, y0] = solve_system(parms);
    yOut{i} = yCyclic;
    % Store the simulation details
    yOut{i}.rho0 = repmat(rho0_vals(i), height(yCyclic),1);
    yOut{i}.id = repmat(i, height(yCyclic),1);
    y0_arr(i,:) = y0;
end

% Concatenate into single table (and reorder with rho0 and t at front)
rho0_ts = vertcat(yOut{:});
rho0_ts = movevars(rho0_ts,{'id','rho0','t'},'before',1);

% Save timeseries to file
folder = '../../data/rho0_sweep';
% Check if the folder exists and, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end
writetable(rho0_ts, fullfile(folder,'timeseries.csv'));

% Save model info to file
defaults = get_default_parms_struct();
defaults = rmfield(defaults,'rho0');
write_sim_info(fullfile(folder,'model_info.yml'),defaults,y0_arr,'rho0_sweep')


%% Heatmap examples
%===============================================================================
% Define the systems of interest
eg_parms = struct( ...
    "A", struct("rho0", 0.1, "omega", 10^-5, "beta1", 10^-6), ...
    "E", struct("rho0", 0.1,  "omega", 10^-5.7, "beta1", 10^-6.0), ...
    "B", struct("rho0", 0.1,  "omega", 10^-5, "beta1", 10^-6.9), ...
    "C", struct("rho0", 0.1,  "omega", 10^-5.25, "beta1", 10^-6.4), ...
    "D", struct("rho0", 0.1,  "omega", 10^-5.18, "beta1", 10^-6.2));
    

% Loop over the system and solve the systems
labels = fields(eg_parms);
resTimeseries = cell(length(labels),1);
y0_arr = zeros(length(labels),8);
for i = 1:length(labels)
    % Get the label and parameters from the structure
    labeli = labels{i}
    parms = eg_parms.(labeli)

    % Solve
    [yCyclic, yConst, y0Cyclic] = solve_system(parms);

    % Add the simulation label to the table
    yCyclic.id = repmat(labeli,height(yCyclic),1);

    % Add to the output structure
    resTimeseries{i} = yCyclic;
    y0_arr(i,:) = y0Cyclic;
end

% Plot all EA
figure;
for i = 1:length(labels)
    subplot(1,length(labels),i);
    labeli = labels{i};
    resi = resTimeseries{i};
    plot(resi.t,resi.EA)
    title(labeli);
end
% Plot all system responses for example B
plot_system_response(resTimeseries{(labels=="E")})


% Save the timeseries data
resTimeseries = vertcat(resTimeseries{:});
resTimeseries = movevars(resTimeseries,{'id','t'},'before',1);
folder = '../../data/examples/';
% Check if the folder exists and, if not, create it
if ~exist(folder, 'dir')
    mkdir(folder);
end
writetable(resTimeseries, fullfile(folder,'timeseries.csv'));

% Save the parameter information
eg_labs = fieldnames(eg_parms);
parms_table = cell(numel(eg_labs),1);
for i = 1:numel(eg_labs)
    tmp = eg_parms.(eg_labs{i});
    tmp.id = eg_labs{i};
    parms_table{i} = struct2table(tmp);
end
defaultparms = get_default_parms_struct();
defaultparms = rmfield(defaultparms,fieldnames(eg_parms.A));
write_sim_info(fullfile(folder,'model_info.yml'), ...
    defaultparms, y0_arr, "heatmap-examples");
writetable(vertcat(parms_table{:}),fullfile(folder,'parameters.csv'))






%% Functions
%===============================================================================

function [yCyclic, yConst, y0Cyclic] = solve_system(parms)
    % Initial conditions
    varNames = {'M0', 'M1',      'M2',  'K0', 'KA','E0','EF','EA'};
    y0       = [9e5;   1e2;     4.5e4;   5e5;  2e4;   0;   0;   0];  
    %y0 = [898042.5; 9637.075; 44814.18; 42929.06; 296376.2; 1260; 105.3528; 107.2753];
    %dY = [-1; 1; -1; -1; -1; 0; 1; 1];
    %y0 = y0.*(1 + 0.1*dY)
    
    % Timing parameters
    t0 = 0;
    dt = 0.5; %0.1;
    tEndConst = 5e3;
    tEndCyclic = 3*13*28;
    
    % DE solver options
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

    % Solve using stiff solver
    [tConst,y1] = ode15s(@(t,y)fn_immune_ode_system(t, y, parms, true), ...
                    (t0:dt:tEndConst), y0', opts);
    y0Cyclic = median(y1);
    [tCyclic,y2] = ode15s(@(t,y)fn_immune_ode_system(t, y, parms, false), ...
                      (t0:dt:tEndCyclic), y0Cyclic, opts);
    
    % Convert the solution to a table and assign column names
    yConst = array2table(y1, 'VariableNames', varNames);
    yConst.t = tConst;

    yCyclic = array2table(y2, 'VariableNames', varNames);
    yCyclic.t = tCyclic;
end


function plot_system_response(y)
    % Easy access to t
    t = y.t;

    % Plot the macrophage (M0, M1, M2)
    figure;
    subplot(1,3,1);
    plot(t, y.M0, '-', 'DisplayName', 'M0');
    hold on;
    plot(t, y.M1, '-', 'DisplayName', 'M1');
    plot(t, y.M2, '-', 'DisplayName', 'M2');
    hold off;
    set(gca, 'YScale', 'log')
    xlabel('Time (days)');
    ylabel('cells/mL');
    title('Macrophages');
    legend show;
    
    % Plot the natural killer cells (K0, KA)
    subplot(1,3,2);
    plot(t, y.K0, '-', 'DisplayName', 'K0');
    hold on;
    plot(t, y.KA, '-', 'DisplayName', 'KA');
    hold off;
    set(gca, 'YScale', 'log')
    xlabel('Time (days)');
    ylabel('cells/mL');
    title('Natural killer cells');
    legend show;
    
    % Plot the endometrial cells (EF, EA)
    subplot(1,3,3);
    %plot(t, y.E0, '-', 'DisplayName', 'E0');
    %hold on;
    plot(t, y.EF, '-', 'DisplayName', 'EF');
    hold on;
    plot(t, y.EA, '-', 'DisplayName', 'EA');
    hold off;
    set(gca, 'YScale', 'log')
    xlabel('Time (days)');
    ylabel('cells/mL');
    title('Endometrial cells');
    legend show;
end

function write_sim_info(filename, parmsMap, icArray, id)
    % File identifier
    fileID = fopen(filename,'w');

    % Write simulation identifier
    fprintf(fileID,"simulation-id: %s\n", id);

    % Write todays date
    fprintf(fileID, 'date: %s\n', datetime('today'));
    
    % Write the default parameter values
    fprintf(fileID,"parameters:\n"); % Title
    parmNames = fieldnames(parmsMap);
    for p=parmNames'
        fprintf( fileID, "\t%s: %.3e\n", p{1}, parmsMap.(p{1}) );
    end

    % Write the initial conditions array (note can be 2 dimensional)
    fprintf(fileID,"initial-conditions:\n"); % Title
    varNames = {'M0', 'M1', 'M2','K0', 'KA','E0','EF','EA'};
    assert(length(varNames)==length(icArray));
    for i=1:length(varNames)
        fprintf(fileID, "\t%s: [ ",varNames{i});
        for j = 1:height(icArray)
            fprintf(fileID, "%.2f\t", icArray(j,i));
        end
        fprintf(fileID, "]\n");
    end
    
    % Close the file
    fclose(fileID);
end

function default_parms_struct = get_default_parms_struct()
    % Generate from script then create structure using table for convenience
    [muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0,deltaE,sigma,MC, ...
     KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK,gamma, ...
        rhoF, muE,omega,beta1,rho0] = default_parameters;
    default_parms_table = table(muM,muK,etaE,etaK,etaM,deltaM,deltaK,deltaE0, ...
        deltaE,sigma,MC,KC,EC,CM1,CM2,CKA,beta2,beta12,beta21,thetaM,thetaK, ...
        gamma, rhoF, muE,omega,beta1,rho0);
    default_parms_struct = table2struct(default_parms_table);
end
