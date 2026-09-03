%% Automated Iterative Sensitivity Analysis and Parameter Optimization


clear all;
close all;
clc;

diaryname = sprintf('_diary-%s-%s.txt', mfilename, datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
diary(diaryname);

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

geometry = '1d';
max_iterations = 3; %
use_equivalent_eff_cond = strcmpi(geometry, '1d');
include_current_collectors = strcmpi(geometry, '1d');
useRegionBruggemanCoefficients = true;
debug = false;
hessian = true;
hessianSteps = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9];

printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));
getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);

%% Configuration

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
    numTimesteps = 400;
    input_ref_1 = struct('I', expdata.I, ...
                         'totalTime', expdata.time(end), ...
                         'numTimesteps', numTimesteps, ...
                         'geometry', geometry, ...
                         'include_current_collectors', include_current_collectors, ...
                         'lowRateParams', jsonEQC, ...
                         'highRateParams', jsonHRC_current_sim, ...
                         'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
    output_ref_1 = runHydra(input_ref_1, 'clearSimulation', false);

    % Evaluate the experimental voltage at the fixed simulation times.
    simtimes = getTime(output_ref_1.states);
    assert(expdata.time(1) <= simtimes(1));
    assert(abs(expdata.time(end) - simtimes(end)) / expdata.time(end) < 1e-14);

    Evals = interp1(expdata.time, expdata.E, simtimes, 'linear', 'extrap');
    statesExp = cell(numel(output_ref_1.states), 1);
    for k = 1:numel(output_ref_1.states)
        statesExp{k}.time = simtimes(k);
        statesExp{k}.(ctrl).E = Evals(k);
    end

    %% Sensitivity Analysis
    fprintf('\n--- Sensitivity Analysis ---\n');

    if useRegionBruggemanCoefficients
        shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', 'elyte_bg_ne', 'elyte_bg_pe', 'elyte_bg_sep'};
    else
        %shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', 'ne_j0', 'pe_j0', 'elyte_bg'};
        shortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', 'elyte_bg'};
    end

    PS = ParamSetter('shortnames', shortnames, ...
                     'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
    if iteration == 0
        PS.printBoxLims();
    end

    simSetup = SimulationSetup(struct('model', output_ref_1.model, ...
                                      'schedule', output_ref_1.schedule, ...
                                      'initstate', output_ref_1.initstate, ...
                                      'NonLinearSolver', output_ref_1.nls, ...
                                      'OutputMinisteps', false));
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

    objFunc = @(simsetup, states, varargin) leastSquaresEI(simsetup, states, statesExp, varargin{:});
    vals = objFunc(simSetup, output_ref_1.states);
    objScaling = sum([vals{:}]);
    objective = @(p, varargin) evalObjectiveBattmo(p, objFunc, simSetup, parameters, ...
                                                   'objScaling', objScaling, varargin{:});

    [v, gradient] = objective(u, 'gradientMethod', 'AdjointAD');

    if debug

        % Compare with finite difference gradient
        pertsize = 1e-6;
        [v_fd, gradient_fd] = objective(u, 'gradientMethod', 'PerturbationADNUM', ...
                                        'PerturbationSize', pertsize);

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
    else
        grouping_strategy = 'hybrid_adaptive';
        % grouping_strategy = 'magnitude';
        % grouping_strategy = 'physical';
        parameter_groups = createClusteredParameterGroups(PS.shortnames, gradient_actual, grouping_strategy, 3);
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
    current_jsonHRC = PS.buildParameterJson(PS.getValues(current_output.model));
    optimization_history = struct();

    for group_idx = 1:length(parameter_groups)
        group = parameter_groups{group_idx};

        fprintf('\nGroup %d: %s\n', group_idx, group.name);

        % Setup optimization
        simSetupGroup = SimulationSetup(struct('model', current_output.model, ...
                                               'schedule', current_output.schedule, ...
                                               'initstate', current_output.initstate, ...
                                               'NonLinearSolver', current_output.nls, ...
                                               'OutputMinisteps', false));

        PS_group = ParamSetter('shortnames', group.parameters, ...
                               'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients);
        PS_group.validate(simSetupGroup.model);

        getValues_group = @(model, notused) PS_group.getValues(model);
        setValues_group = @(model, notused, v) PS_group.setValues(model, v);

        params_group{1} = ModelParameterScaled(simSetupGroup, ...
                                               'name', 'ParamSetter', ...
                                               'belongsTo', 'model', ...
                                               'location', {''}, ...
                                               'boxLims', PS_group.boxLims, ...
                                               'scaling', 'nan', ...
                                               'shortnames', PS_group.shortnames, ...
                                               'getfun', getValues_group, ...
                                               'setfun', setValues_group);

        X0 = getScaledParameterVector(simSetupGroup, params_group);

        % Objective function with scaling
        objFuncGroup = @(simsetup, states, varargin) leastSquaresEI(simsetup, states, statesExp, varargin{:});
        vals_group = objFuncGroup(simSetupGroup, current_output.states);

        objScaling_group = sum([vals_group{:}]);

        % Optimization
        objectiveGradient = @(p, varargin) evalObjectiveBattmo(p, objFuncGroup, simSetupGroup, params_group, ...
                                                               'objScaling', objScaling_group, varargin{:});

        gradTol = 1e-3;
        objChangeTol = -inf;
        maxit = 500;
        [vOpt, Xopt, history] = unitBoxBFGS(X0, objectiveGradient, ...
                                            'gradTol'         , gradTol     , ...
                                            'objChangeTol'    , objChangeTol, ...
                                            'lineSearchMaxIt' , 10          , ...
                                            'maxInitialUpdate', 0.02        , ...
                                            'maximize'        , false       , ...
                                            'maxit'           , maxit       , ...
                                            'logPlot'         , true        , ...
                                            'limitedMemory'   , ~hessian    , ...
                                            'outputHessian'   , hessian);
        reasonStr = getReasonStr(history, ...
                                 'gradTol'     , gradTol     , ...
                                 'objChangeTol', objChangeTol, ...
                                 'maxit'       , maxit);
        disp(reasonStr);
        disp(getBoxLimHits(simSetupGroup, params_group, PS_group, Xopt));

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
        setupOpt = updateSetupFromScaledParameters(simSetupGroup, params_group, Xopt);

        current_jsonHRC = PS.buildParameterJson(PS.getValues(setupOpt.model));

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
    disp('HRC params:');
    printer(current_jsonHRC);
    fprintf('End of iteration %g\n', iteration);

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

if any(strcmp(PS_final.shortnames, 'elyte_bgfactor'))
    final_model = current_output.model;
    optimal_bgfactor = final_model.(elyte).bgfactor;
    region_bruggeman = final_model.(elyte).regionBruggemanCoefficients;

    effective_region_bruggeman = struct( ...
        ne, optimal_bgfactor * region_bruggeman.(ne), ...
        pe, optimal_bgfactor * region_bruggeman.(pe), ...
        sep, optimal_bgfactor * region_bruggeman.(sep));

    fprintf('\nOptimal electrolyte Bruggeman factor: %.16g\n', optimal_bgfactor);
    disp('Effective electrolyte region Bruggeman coefficients:');
    printer(effective_region_bruggeman); % Can be the same value if the same bg initial guess is used

    % Tortuosities
    tags = final_model.(elyte).regionTags;
    poro = struct();
    poro.(ne) = unique(final_model.(elyte).volumeFraction(tags == 1));
    poro.(pe) = unique(final_model.(elyte).volumeFraction(tags == 2));
    poro.(sep) = unique(final_model.(elyte).volumeFraction(tags == 3));
    assert(numel(poro.(ne)) == 1);
    assert(numel(poro.(pe)) == 1);
    assert(numel(poro.(sep)) == 1);

    effective_bruggeman_json = struct();
    effective_bruggeman_json.(elyte).regionBruggemanCoefficients = effective_region_bruggeman;
    final_tortuosities = calculateTortuosityFromBruggeman(poro, effective_bruggeman_json);

    disp('Final electrolyte tortuosities:');
    printer(final_tortuosities);
end

%% Check hessian for the final optimization group

if hessian
    invHscaled = history.hess{end};
    hessianfdpertsize = 1e-6;

    debug = true;
    Hscaled = calculateBFGSHessian(invHscaled, PS_group.shortnames);
    if debug
        [HfdComparison, HfdReport] = calculateFDHessian(Xopt, objectiveGradient, PS_group.shortnames, hessianSteps);
        compareHessians(Hscaled, HfdComparison, HfdReport, Xopt, PS_group.shortnames);
    end
    HfdScaled = calculateFDHessian(Xopt, objectiveGradient, PS_group.shortnames, hessianfdpertsize);

    plotHessianEigenvectors(Hscaled, PS_group.shortnames, 'BFGS');
    plotHessianEigenvectors(HfdScaled, PS_group.shortnames, 'FD');
end

diary off;
