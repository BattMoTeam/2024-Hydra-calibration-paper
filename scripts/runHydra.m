function output = runHydra(input, varargin)

    % Input parameters
    input_default = struct('DRate'                         , []   , ...
                           'totalTime'                     , []   , ...
                           'numTimesteps'                  , 100  , ...
                           'lowRateParams'                 , []   , ...
                           'highRateParams'                , []   , ...
                           'useExpDiffusion'               , false, ...
                           'neD'                           , []   , ...
                           'peD'                           , []   , ...
                           'useRegionBruggemanCoefficients', false, ...
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
                 'validateJson'   , true);

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

    % Load base json
    jsonstruct = parseBattmoJson(fullfile(getHydra0Dir(), 'parameters', 'h0b-base.json'));

    jsonstruct.include_current_collectors = input.include_current_collectors;

    % Use experimental diffusion if requested
    if input.useExpDiffusion
        assert(isempty(input.highRateParams), ...
               'Should not use experimental diffusion and high rate params simultaneously');

        eldes = {ne, pe};
        for ielde = 1:numel(eldes)
            elde = eldes{ielde};

            % Remove existing diffusion params
            assert(isfield(jsonstruct.(elde).(co).(am).(sd), 'referenceDiffusionCoefficient'), ...
                   'Expected diffusion parameters to be present in the base json');
            jsonstruct.(elde).(co).(am).(sd) = rmfield(jsonstruct.(elde).(co).(am).(sd), 'referenceDiffusionCoefficient');

            % Add experimental diffusion params
            switch elde
              case ne
                functionname = 'computeDanodeH0b';
              case pe
                functionname = 'computeDcathodeH0b';
              otherwise
                error('Unexpected electrode %s', elde);
            end
            jsonstruct_diffusion = struct('type', 'function', ...
                                          'functionname', functionname, ...
                                          'argumentlist', 'soc');
            jsonstruct.(elde).(co).(am).(sd).diffusionCoefficient = jsonstruct_diffusion;
        end
    end

    % Set input diffusion
    if not(isempty(input.neD))
        jsonstruct.(ne).(co).(am).(sd).referenceDiffusionCoefficient = input.neD;
    end

    if not(isempty(input.peD))
        jsonstruct.(pe).(co).(am).(sd).referenceDiffusionCoefficient = input.peD;
    end

    % Set low rate params
    if not(isempty(input.lowRateParams))
        jsonstruct_low_rate_params = input.lowRateParams;
        jsonstruct = mergeJsonStructs({jsonstruct_low_rate_params, jsonstruct}, 'warn', false);
    end

    % Set high rate params
    if not(isempty(input.highRateParams))
        jsonstruct_high_rate_params = input.highRateParams;
        jsonstruct = mergeJsonStructs({jsonstruct_high_rate_params, jsonstruct}, 'warn', false);
    end

    if input.useRegionBruggemanCoefficients
        jsonstruct.(elyte).useRegionBruggemanCoefficients = true;
        jsonstruct.(elyte).regionBruggemanCoefficients = struct(ne, 1.5, ...
                                                                pe, 1.5, ...
                                                                sep, 1.5);
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
    jsonstruct = mergeJsonStructs({jsonstruct_geom, jsonstruct});

    % Scale input geometry
    if strcmpi(jsonstruct.Geometry.case, "1D") && jsonstruct.include_current_collectors
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
    model = Battery(paramobj);

    % Setup nonlinear solver
    jsonstruct_nls = parseBattmoJson(fullfile('Utilities', 'Linearsolvers', 'JsonDataFiles', 'default_direct_linear_solver.json'));
    jsonstruct_nls.verbose = opt.verbose;
    jsonstruct = mergeJsonStructs({jsonstruct_nls, jsonstruct});
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

    dt = totalTime / input.numTimesteps;
    dt = rampupTimesteps(totalTime, dt, 10, 'threshold_error', 1e-8);
    step = struct('val', dt, 'control', ones(numel(dt), 1));

    tup = 1*minute;
    cutOffVoltage = model.(ctrl).lowerCutoffVoltage;
    srcfunc = @(t, I, E, Imax) rampupSwitchControl(t, tup, I, E, model.(ctrl).Imax, cutOffVoltage);
    control = struct('src', srcfunc);
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
                                               'NonLinearSolver', nls);
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
            [~, output.states] = getPackedSimulatorOutput(output.problem);

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
        writeJsonStruct(jsonencode(input, 'PrettyPrint', true), jsoninputfilename);

        if opt.clearSimulation
            clearPackedSimulatorOutput(output.problem, 'Prompt', false);
        end

        simulatePackedProblem(output.problem);

        if nargout > 0
            [~, output.states] = getPackedSimulatorOutput(output.problem);
        end

    else

        [~, output.states] = simulateScheduleAD(initstate, model, schedule, ...
                                                'OutputMinisteps', true, ...
                                                'NonLinearSolver', nls);

    end

    output.input = input;

end

%{
  Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
