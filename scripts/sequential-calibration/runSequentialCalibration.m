%% Automated Iterative Sensitivity Analysis and Parameter Optimization


clear all;
close all;
clc;

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
itf     = 'Interface';
elyte   = 'Electrolyte';
sep     = 'Separator';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
cc      = 'CurrentCollector';
ctrl    = 'Control';
sd      = 'SolidDiffusion';
co      = 'Coating';
bd      = 'Binder';
ad      = 'ConductingAdditive';

diary(sprintf('diary-runSquentialCalibration-%s.txt', datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-ss')));

%% Configuration
geometry = '1d';
max_iterations = 3;
use_equivalent_eff_cond = strcmpi(geometry, '1d');
include_current_collectors = strcmpi(geometry, '1d');
useRegionBruggemanCoefficients = false;

% Experimental data and base parameters
datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Lowest DRate is first
expdata_all = cell(1, numel(dataraw.time));
for idata = 1:numel(dataraw.time)
    time = dataraw.time{idata} * hour;
    current = dataraw.current{idata};
    expdata_all{idata} = struct('time', time                         , ...
                                'U'   , dataraw.voltage{idata}       , ...
                                'E'   , dataraw.voltage{idata}       , ...
                                'cap' , abs(trapz(time, current))    , ...
                                'I'   , abs(mean(current))           , ...
                                'rawI', current                      , ...
                                'Q'   , abs(cumtrapz(time, current)));
    expdata_all{idata}.DRate = expdata_all{idata}.I / expdata_all{idata}.cap * hour;
    expdata_all{idata}.rate = expdata_all{idata}.DRate;
end

% Highest DRate is last
expdata = expdata_all{end};
optimization_lower_cutoff_voltage = min(expdata.E) - 0.5;
jsonEQC = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json'));

% Utility functions
getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);


%% Main Optimization Loop
fprintf('=== ITERATIVE PARAMETER OPTIMIZATION ===\n');

% Store results for all iterations
all_iteration_results = struct();
initial_calibrated_values = [];
initial_group_by_shortname = {};

for iteration = 0:max_iterations
    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('ITERATION %d\n', iteration);
    fprintf('%s\n', repmat('=', 1, 60));

    %% Load parameters
    if iteration == 0
        jsonHRC_current = [];
    else
        jsonHRC_current = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', ...
                                                   sprintf('h0b-%s-iterative-%d-groups-%s.json', grouping_strategy, iteration, geometry)));
    end

    jsonHRC_current_sim = jsonHRC_current;
    if ~isempty(jsonHRC_current_sim)
        jsonHRC_current_sim.(ctrl).lowerCutoffVoltage = optimization_lower_cutoff_voltage;
    end

    %% Run simulation
    input_ref_1 = struct('I', expdata.I, ...
                         'totalTime', expdata.time(end), ...
                         'geometry', geometry, ...
                         'include_current_collectors', include_current_collectors, ...
                         'lowRateParams', jsonEQC, ...
                         'highRateParams', jsonHRC_current_sim, ...
                         'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
    output_ref_1 = runHydra(input_ref_1, 'clearSimulation', false);

    % % Truncate if needed
    % if output_ref_1.states{end}.time > expdata.time(end)
    %     idx = find(getTime(output_ref_1.states) <= expdata.time(end), 1, 'last');
    %     output_ref_1.states = output_ref_1.states(1:idx);
    % end

    %% Sensitivity Analysis
    fprintf('\n--- Sensitivity Analysis ---\n');

    if useRegionBruggemanCoefficients
        shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', ...%'ne_j0', 'pe_j0', ...
                      'elyte_bg_ne', 'elyte_bg_pe', 'elyte_bg_sep'};
        % shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', 'elyte_bg_ne', 'elyte_bg_pe', 'elyte_bg_sep'};
    else
        %shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', 'ne_j0', 'pe_j0', 'elyte_bg'};
        shortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', 'elyte_bg'};
    end

    PS = ParamSetter('shortnames', shortnames, ...
                     'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
    if iteration == 0
        PS.printBoxLims();
    end

    simSetup = struct('state0', output_ref_1.initstate, ...
                      'model', output_ref_1.model, ...
                      'schedule', output_ref_1.schedule);
    PS.validate(simSetup.model);

    if iteration == 0
        initial_calibrated_values = PS.getValues(simSetup.model);
    end

    getValues = @(model, notused) PS.getValues(model);
    setValues = @(model, notused, v) PS.setValues(model, v);

    parameters{1} = ModelParameterScaled(simSetup, ...
                                         'name', 'ParamSetter', ...
                                         'belongsTo', 'model', ...
                                         'location', {''}, ...
                                         'boxLims', PS.boxLims, ...
                                         'scaling', 'linear', ...
                                         'shortnames', PS.shortnames, ...
                                         'getfun', getValues, ...
                                         'setfun', setValues);

    % Objective function with proper scaling
    u = getScaledParameterVector(simSetup, parameters);

    warning('check paramoptfuncreg objective');
    paramoptfunc = ParamOptFuncReg(expdata, PS.shortnames, simSetup, u);
    objFunc = @(model, states, schedule, parameters,p, varargin) ...
              paramoptfunc.leastSquaresEWithReg(model, states, schedule, parameters,p, varargin{:});

    vals = objFunc(simSetup.model, output_ref_1.states, simSetup.schedule, parameters, u);
    objScaling = sum([vals{:}]);

    [v, gradient, ~] = evalObjectiveBattmoReg(...
        u, objFunc, simSetup, parameters, ...
        'gradientMethod', 'AdjointAD', 'objScaling', objScaling);

    debug = false;
    if debug

        % Compare with finite difference gradient
        pertsize = 1e-6;
        [v_fd, gradient_fd] = evalObjectiveBattmoReg(...
            u, objFunc, simSetup, parameters, ...
            'gradientMethod', 'PerturbationADNUM', ...
            'PerturbationSize', pertsize, ...
            'objScaling', objScaling);

        relErr = abs(gradient - gradient_fd) ./ (abs(gradient) + abs(gradient_fd));
        tol = 1e-4;
        isLarge = char((relErr > tol) * '*');
        disp(table(PS.shortnames(), gradient, gradient_fd, relErr, isLarge));
        assert(abs(v - v_fd) < eps);
        assert(all(abs(gradient) > 0));
        assert(all(abs(gradient_fd) > 0));
        assert(norm(relErr, 'inf') < tol);

        return
    end

    gradient_actual = gradient * objScaling;

    if isempty(gradient_actual)
        fprintf('Skipping iteration %d (empty gradient)\n', iteration);
        continue;
    end

    fprintf('Parameter Sensitivities:\n');
    for i = 1:length(PS.shortnames)
        fprintf('  %-10s: %.3e\n', PS.shortnames{i}, abs(gradient_actual(i)));
    end

    %% Parameter Grouping
    fprintf('\n--- Parameter Grouping ---\n');
    parameter_groups = {};

    % Alternate between eg hybrid adaptive grouping strategy and
    % optimize all parameters
    if mod(iteration, 2)
        parameter_groups{1} = struct(...
            'name', 'All_Parameters', ...
            'parameters', {PS.shortnames}, ...
            'priority', 1);
        OptimizationSolver = 'unitboxbfgs';
    else
        grouping_strategy = 'hybrid_adaptive';
        %grouping_strategy = 'magnitude'; % cannot run validation: conv problems
        parameter_groups = createClusteredParameterGroups(PS.shortnames, gradient_actual, grouping_strategy, 3);
        OptimizationSolver = 'unitboxbfgs';
    end

    if iteration == 0
        initial_group_by_shortname = repmat({''}, numel(PS.shortnames), 1);
        for group_index = 1:numel(parameter_groups)
            group = parameter_groups{group_index};
            for parameter_index = 1:numel(group.parameters)
                match = strcmp(PS.shortnames, group.parameters{parameter_index});
                assert(any(match), 'Group parameter %s is not in the complete shortname list.', group.parameters{parameter_index});
                initial_group_by_shortname{match} = group.name;
            end
        end
    end

    fprintf('Groups for optimization:\n');
    for i = 1:length(parameter_groups)
        fprintf('  %d. %s: %s\n', i, parameter_groups{i}.name, strjoin(parameter_groups{i}.parameters, ', '));
    end

    %% Sequential Group Optimization
    fprintf('\n--- Group Optimization ---\n');

    current_output = output_ref_1;
    current_jsonHRC = buildParameterJson(PS, PS.getValues(current_output.model));
    optimization_history = struct();

    for group_idx = 1:length(parameter_groups)
        group = parameter_groups{group_idx};

        fprintf('\nGroup %d: %s\n', group_idx, group.name);

        % Store state for plotting
        output_before = current_output;

        % Setup optimization
        simSetup_group = struct('state0', current_output.initstate, ...
                                'model', current_output.model, ...
                                'schedule', current_output.schedule);

        PS_group = ParamSetter('shortnames', group.parameters, ...
                               'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
        PS_group.validate(simSetup_group.model);

        getValues_group = @(model, notused) PS_group.getValues(model);
        setValues_group = @(model, notused, v) PS_group.setValues(model, v);

        params_group{1} = ModelParameterScaled(simSetup_group, ...
                                               'name', 'ParamSetter', ...
                                               'belongsTo', 'model', ...
                                               'location', {''}, ...
                                               'boxLims', PS_group.boxLims, ...
                                               'scaling', 'nan', ...
                                               'shortnames', PS_group.shortnames, ...
                                               'getfun', getValues_group, ...
                                               'setfun', setValues_group);

        X0 = getScaledParameterVector(simSetup_group, params_group);

        % Objective function with scaling
        paramoptfunc_group = ParamOptFuncReg(expdata, group.parameters, simSetup_group, X0);
        objFunc_group = @(model, states, schedule, params_group, varargin) ...
            paramoptfunc_group.leastSquaresEWithReg(model, states, schedule, params_group, varargin{:});
        vals_group = objFunc_group(simSetup_group.model, current_output.states, simSetup_group.schedule, params_group,X0);

        objScaling_group = sum([vals_group{:}]);

        % Optimization
        objectiveGradient = @(p) evalObjectiveBattmoReg(p, objFunc_group, simSetup_group, params_group, ...
                                                        'objScaling', objScaling_group, 'OptimizationSolver', OptimizationSolver);

        BFGSopts = {'objChangeTol', 1e-10, 'gradTol', 1e-10, 'maxit', 100, 'maximize', false, 'logPlot', true};
        [vOpt, Xopt, history] = unitBoxBFGS(X0, objectiveGradient, BFGSopts{:});

        disp(getBFGSstopReason(history, BFGSopts));
        disp(getBoxLimHits(simSetup_group, params_group, PS_group, Xopt));

        % Store results
        optimization_history(group_idx).group_name = group.name;
        optimization_history(group_idx).parameters = group.parameters;
        optimization_history(group_idx).initial_objective = history.val(1);
        optimization_history(group_idx).final_objective = vOpt;
        optimization_history(group_idx).improvement = (history.val(1) - vOpt) / history.val(1) * 100;
        optimization_history(group_idx).iteration_history = history.val;

        fprintf('  Improvement: %.1f%% (%.3e -> %.3e)\n', ...
                optimization_history(group_idx).improvement, history.val(1), vOpt);
        fprintf('When used group %d: %s\n', group_idx, group.name);
        disp(PS_group.shortnames');

        % Apply optimized parameters
        setupOpt = updateSetupFromScaledParameters(simSetup_group, params_group, Xopt);

        current_jsonHRC = buildParameterJson(PS, PS.getValues(setupOpt.model));

        % Run simulation with new parameters
        current_jsonHRC_sim = current_jsonHRC;
        current_jsonHRC_sim.(ctrl).lowerCutoffVoltage = optimization_lower_cutoff_voltage;
        input_ref_1.highRateParams = current_jsonHRC_sim;
        current_output = runHydra(input_ref_1, 'clearSimulation', false);
    end

    %% Save and Evaluate
    fprintf('\n--- Saving Results ---\n');

    output_filename = sprintf('h0b-%s-iterative-%d-groups-%s.json', grouping_strategy, iteration + 1, geometry);
    writeStruct(current_jsonHRC, fullfile(getHydra0Dir(), 'parameters', output_filename));

    initial_wL2 = l2error(getTime(output_ref_1.states), getE(output_ref_1.states), ...
                          expdata.time, expdata.E, 'extrap', true);
    final_wL2 = l2error(getTime(current_output.states), getE(current_output.states), ...
                        expdata.time, expdata.E, 'extrap', true);

    fprintf('RMSE: %2.2f -> %2.2f mV (%.1f%% reduction)\n', initial_wL2/milli, final_wL2/milli, (initial_wL2-final_wL2)/initial_wL2*100);

    % Store iteration results
    all_iteration_results(iteration+1).iteration = iteration;
    all_iteration_results(iteration+1).initial_wL2 = initial_wL2;
    all_iteration_results(iteration+1).final_wL2 = final_wL2;
    all_iteration_results(iteration+1).improvement = (initial_wL2-final_wL2)/initial_wL2*100;
    all_iteration_results(iteration+1).jsonHRC = current_jsonHRC;

    %% Plot results
    plotIterationResults(iteration, output_ref_1, current_output, expdata, gradient_actual, PS.shortnames, optimization_history);

    %% Validation
    fprintf('\n--- Multi-Rate Validation ---\n');
    plotMultiRateValidation(expdata_all, jsonEQC, current_jsonHRC, iteration, ...
                            'geometry', geometry, ...
                            'include_current_collectors', include_current_collectors, ...
                            'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);

    % %% Check convergence
    % if iteration > 0 && abs((initial_wL2 - final_wL2) / initial_wL2) < 0.01
    %     fprintf('Convergence reached (<1%% improvement)\n');
    %     break;
    % end

end

%% Final calibrated parameter values
final_jsonHRC = all_iteration_results(end).jsonHRC;
PS_final = ParamSetter('shortnames', shortnames, ...
                       'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
final_json_sim = mergeStructs({final_jsonHRC, jsonEQC});
final_calibrated_values = PS_final.getValues(final_json_sim);
fprintf('%.16g\n', final_calibrated_values);

disp(table(PS_final.shortnames(), initial_calibrated_values, final_calibrated_values, ...
           initial_group_by_shortname, ...
           'VariableNames', {'Shortname', 'InitialValue', 'FinalValue', 'InitialGroup'}))

diary off;

function jsonstruct = buildParameterJson(parameterSetter, values)
    locs = parameterSetter.locations();
    jsonstruct = struct();
    for index = 1:numel(locs)
        jsonstruct = setStructField(jsonstruct, locs{index}, values(index));
    end
end
