function output = runHydra(input, varargin)

    % Input parameters
    input_default = struct('I'                             , []   , ...
                           'DRate'                         , []   , ...
                           'totalTime'                     , []   , ...
                           'numTimesteps'                  , 100  , ...
                           'lowRateParams'                 , []   , ...
                           'highRateParams'                , []   , ...
                           'useRegionBruggemanCoefficients', false, ...
                           'OutputMinisteps'               , false, ...
                           'include_current_collectors'    , false, ...
                           'geometry'                      , '1d');

    if not(isempty(input))
        fds = fieldnames(input);
        vals = cellfun(@(fd) input.(fd), fds, 'un', false);
        input = horzcat(fds, vals);
        input = reshape(input', [], 1);
        input = merge_options(input_default, input{:});
    else
        input = input_default;
    end

    % Solver options
    opt = struct('runSimulation'  , true    , ...
                 'dopacked'       , true    , ...
                 'verbose'        , false   , ...
                 'clearSimulation', true    , ...
                 'outputDirectory', 'output', ...
                 'validateJson'   , false);

    opt = merge_options(opt, varargin{:});

    % Handy short names
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    elyte = 'Electrolyte';
    am    = 'ActiveMaterial';
    itf   = 'Interface';
    sd    = 'SolidDiffusion';
    ctrl  = 'Control';
    co    = 'Coating';
    bd    = 'Binder';
    ca    = 'ConductingAdditive';
    cc    = 'CurrentCollector';
    sep   = 'Separator';
    geom  = 'Geometry';
    rbc   = 'regionBruggemanCoefficients';

    % Load base json
    jsonstruct = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));
    jsonstruct.include_current_collectors = input.include_current_collectors;

    % Set low rate params
    if not(isempty(input.lowRateParams))
        jsonstruct_low_rate_params = input.lowRateParams;
        jsonstruct = mergeStructs({jsonstruct_low_rate_params, jsonstruct}, 'warn', false);
    end

    % Set high rate params
    if not(isempty(input.highRateParams))
        jsonstruct_high_rate_params = input.highRateParams;
        jsonstruct = mergeStructs({jsonstruct_high_rate_params, jsonstruct}, 'warn', false);
    end

    if input.useRegionBruggemanCoefficients
        jsonstruct.(elyte).useRegionBruggemanCoefficients = true;

        bgFromTau = @(poro, tau) -log(tau) / log(poro);

        % Set if not already set (via jsonstructHRC)
        if ~isfield(jsonstruct.(elyte), rbc)
            jsonstruct.(elyte).regionBruggemanCoefficients = struct();
        end
        if ~isfield(jsonstruct.(elyte).(rbc), ne)
            poro = 1 - jsonstruct.(ne).(co).volumeFraction;
            tauref = 3.46;
            jsonstruct.(elyte).regionBruggemanCoefficients.(ne) = bgFromTau(poro, tauref);
        end
        if ~isfield(jsonstruct.(elyte).(rbc), pe)
            poro = 1 - jsonstruct.(pe).(co).volumeFraction;
            tauref = 3.;
            jsonstruct.(elyte).regionBruggemanCoefficients.(pe) = bgFromTau(poro, tauref);
        end
        if ~isfield(jsonstruct.(elyte).(rbc), sep)
            poro = jsonstruct.(sep).porosity;
            tauref = 4.2;
            jsonstruct.(elyte).regionBruggemanCoefficients.(sep) = bgFromTau(poro, tauref);
        end
    end

    % Load geometry
    switch lower(input.geometry)
      case '1d'
        geomfile = 'h0b-geometry-1d.json';
      case '3d'
        geomfile = 'h0b-geometry-3d.json';
      otherwise
        error('Unsupported geometry %s', input.geometry);
    end
    jsonstruct_geom = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', geomfile));
    jsonstruct = mergeStructs({jsonstruct_geom, jsonstruct});

    % Scale input geometry
    if strcmpi(jsonstruct.Geometry.case, '1D') && jsonstruct.include_current_collectors
        json_geom_3d = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-geometry-3d.json'));

        ne_LH = struct('L', json_geom_3d.(geom).length, ...
                       'h', json_geom_3d.(geom).width, ...
                       't', jsonstruct.(ne).(cc).thickness);

        ne_effkappa = geometryScaling(ne_LH, jsonstruct.(ne).(cc).electronicConductivity);
        jsonstruct.(ne).(cc).electronicConductivity = ne_effkappa;

        pe_LH = struct('L', json_geom_3d.(geom).length, ...
                       'h', json_geom_3d.(geom).width, ...
                       't', jsonstruct.(pe).(cc).thickness);

        pe_effkappa = geometryScaling(pe_LH, jsonstruct.(pe).(cc).electronicConductivity);
        jsonstruct.(pe).(cc).electronicConductivity = pe_effkappa;

    end

    % Validate json (requires python)
    if opt.validateJson
        validateJsonStruct(jsonstruct);
    end

    % Convert to battery input parameters
    paramobj = BatteryInputParams(jsonstruct);
    paramobj = setupBatteryGridFromJson(paramobj, jsonstruct);

    % Set rate if provided
    if not(isempty(input.DRate))
        paramobj.(ctrl).DRate = input.DRate;
    end

    % Validate before building model
    paramobj = paramobj.validateInputParams();
    model = GenericBattery(paramobj);

    % Set rate or current
    hasDRate = not(isempty(input.DRate));
    hasI = not(isempty(input.I));
    if hasDRate && hasI
        error('Specify either DRate or I, not both');
    end
    if not(hasDRate) && not(hasI) && opt.runSimulation
        error('Specify either DRate or I');
    end
    if hasI
        model.(ctrl).Imax = input.I;
    end
    if hasDRate
        model.(ctrl).DRate = input.DRate;
    end

    % Setup nonlinear solver
    jsonstruct_nls = parseBattmoJson(fullfile('Utilities', 'Linearsolvers', 'JsonDataFiles', 'default_direct_linear_solver.json'));
    jsonstruct_nls.verbose = opt.verbose;
    jsonstruct = mergeStructs({jsonstruct_nls, jsonstruct});
    [model, nls, jsonstruct] = setupNonLinearSolverFromJson(model, jsonstruct);

    % Basic config
    model.verbose = opt.verbose;
    model.AutoDiffBackend = AutoDiffBackend();

    % Setup initial state and time stepping
    initstate = model.setupInitialState();

    if isempty(input.totalTime)
        totalTime = 1*hour / model.(ctrl).DRate;
    else
        totalTime = input.totalTime;
    end

    timestep = struct('totalTime', totalTime, ...
                      'numberOfTimeSteps', input.numTimesteps, ...
                      'useRampup', true, ...
                      'numberOfRampupSteps', 10);
    step    = model.Control.setupScheduleStep(timestep);
    control = model.Control.setupScheduleControl();
    schedule = struct('control', control, 'step', step);

    % Store variables
    output.model      = model;
    output.schedule   = schedule;
    output.paramobj   = paramobj;
    output.initstate  = initstate;
    output.nls        = nls;
    output.jsonstruct = jsonstruct;

    % Setup simulation
    if opt.dopacked
        input.simtag = md5sum(input);

        directory = fullfile(getHydra0Dir(), opt.outputDirectory);
        dataFolder = input.simtag;
        output.problem = packSimulationProblem(initstate, model, schedule, dataFolder, ...
                                               'Directory', directory                , ...
                                               'Name', input.simtag                  , ...
                                               'NonLinearSolver', nls                , ...
                                               'ExtraArguments', ...
                                               {'OutputMinisteps', input.OutputMinisteps});
        output.dataDirectory = output.problem.OutputHandlers.states.dataDirectory;
        output.dataFolder    = output.problem.OutputHandlers.states.dataFolder;
        inputfilename        = fullfile(output.dataDirectory, output.dataFolder, 'input.mat');
        jsoninputfilename    = fullfile(output.dataDirectory, output.dataFolder, 'input.json');

        if not(isempty(input.lowRateParams))
            output.jsonstruct_low_rate_params = jsonstruct_low_rate_params;
        end

        if not(isempty(input.highRateParams))
            output.jsonstruct_high_rate_params = jsonstruct_high_rate_params;
        end

        if not(opt.runSimulation)

            output.input = input;
            [~, output.states, output.reports] = getPackedSimulatorOutput(output.problem);

            if isempty(output.states)
                foundresults = false;
            else
                foundresults = true;
            end

            if foundresults
                dispif(opt.verbose, sprintf('Results of a previous simulation have been found and added to the output\n'));
            elseif opt.verbose
                fprintf('No previous simulations with hash %s were found for this setup in the %s directory\n', ...
                        input.simtag, opt.outputDirectory);
            end

            return

        end

        save(inputfilename, 'input');
        writeStruct(jsonencode(input, 'PrettyPrint', true), jsoninputfilename);

        if opt.clearSimulation
            clearPackedSimulatorOutput(output.problem, 'Prompt', false);
        end

        simulatePackedProblem(output.problem);

        if nargout > 0
            [~, output.states, output.reports] = getPackedSimulatorOutput(output.problem);
        end

    else

        [~, output.states, output.reports] = simulateScheduleAD( ...
            initstate, model, schedule, ...
            'OutputMinisteps', input.OutputMinisteps, ...
            'NonLinearSolver', nls);

    end

    output.input = input;

end

%{
  Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
  and SINTEF Digital, Mathematics & Cybernetics.

  This file is part of The Battery Modeling Toolbox BattMo

  BattMo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BattMo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
