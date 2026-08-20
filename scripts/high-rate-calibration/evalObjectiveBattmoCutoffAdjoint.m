function [objValue, scaledGradient, states, setupNew] = ...
    evalObjectiveBattmoCutoffAdjoint(pvec, objFunc, setup, parameters, varargin)
%EVALOBJECTIVEBATTMOCUTOFFADJOINT Adjoint objective for cutoff simulations.
%
% Run the forward problem with its full schedule so the stop function can
% terminate it at the lower voltage cutoff.  The unused part of the
% schedule is then removed before the objective and reverse solve are
% evaluated.  The crossing state itself is also discarded, so the retained
% domain contains only states strictly above cutoff.  This makes the state
% and schedule lengths consistent for the existing BattMo adjoint
% implementation.
%
% The resulting derivative is conditional on the current terminal time
% step.  It does not include a derivative of the discrete cutoff-step index.
% If the retained states do not reach the required minimum simulation time,
% no admissible objective domain exists.  Return a failed objective
% evaluation so the optimization line search can shorten its trial step and
% try again.

    opt = struct('AdjointLinearSolver', [], ...
                 'objScaling'         , 1 , ...
                 'enforceBounds'      , true, ...
                 'minimumSimulationTime', 0);
    opt = merge_options(opt, varargin{:});

    validateattributes(opt.minimumSimulationTime, {'numeric'}, ...
                       {'scalar', 'real', 'finite', 'nonnegative'});

    pvec = pvec(:);
    numberOfParameters = cellfun(@(parameter) parameter.nParam, parameters);

    if opt.enforceBounds
        pvec = max(0, min(1, pvec));
    end

    parameterVectors = mat2cell(pvec, numberOfParameters, 1);
    parameterValues = cell(size(parameters));
    setupNew = setup;

    for index = 1:numel(parameters)
        parameterValues{index} = parameters{index}.unscale(parameterVectors{index});
        setupNew = parameters{index}.setParameter(setupNew, parameterValues{index});
    end

    states = setupNew.run();
    states = truncateStatesAtCutoff( ...
        states, setupNew.model.Control.lowerCutoffVoltage);

    reachesMinimumTime = false;
    if ~isempty(states)
        finalSimulationTime = states{end}.time;
        timeTolerance = 128 * eps(max([1, abs(finalSimulationTime), ...
                                      abs(opt.minimumSimulationTime)]));
        reachesMinimumTime = finalSimulationTime + timeTolerance >= ...
                             opt.minimumSimulationTime;
    end

    if ~reachesMinimumTime
        objValue = NaN;
        if nargout >= 2
            scaledGradient = nan(size(pvec));
        end
        return
    end

    setupNew = truncateSetupSchedule(setupNew, states);

    objectiveValues = objFunc(setupNew, states);
    assert(all(isfinite([objectiveValues{:}])), ...
           'Objective function values are not finite');
    objValue = sum(vertcat(objectiveValues{:})) / opt.objScaling;
    assert(isfinite(objValue), 'Objective function value is not finite');

    if nargout < 2
        return
    end

    parameterNames = applyFunction(@(parameter) parameter.name, parameters);
    getObjective = @(timeStep, model, state, computeStatePartial) ...
        objectiveAtTimeStep(objFunc, setupNew, states, timeStep, model, ...
                            state, computeStatePartial);

    sensitivities = computeSensitivitiesAdjointADBattmo( ...
        setupNew, states, parameters, getObjective, ...
        'LinearSolver', opt.AdjointLinearSolver);

    gradient = cellfun(@(name) sensitivities.(name), parameterNames, ...
                       'UniformOutput', false);
    scaledGradientParts = cell(numel(parameterNames), 1);

    for index = 1:numel(parameterNames)
        scaledGradientParts{index} = parameters{index}.scaleGradient( ...
            gradient{index}, parameterValues{index});
    end

    scaledGradient = vertcat(scaledGradientParts{:}) / opt.objScaling;

end

function setup = truncateSetupSchedule(setup, states)

    numberOfStates = numel(states);
    numberOfSteps = numel(setup.schedule.step.val);

    assert(numberOfStates > 0, 'The simulation did not produce any states');
    assert(numberOfStates <= numberOfSteps, ...
           'There are more simulated states than scheduled time steps');

    stepFields = fieldnames(setup.schedule.step);
    for index = 1:numel(stepFields)
        field = stepFields{index};
        value = setup.schedule.step.(field);
        if isvector(value) && numel(value) == numberOfSteps
            setup.schedule.step.(field) = value(1:numberOfStates);
        end
    end

    scheduleTime = cumsum(setup.schedule.step.val(:));
    stateTime = cellfun(@(state) state.time, states(:));
    tolerance = 1e-10 * max(1, scheduleTime(end));
    assert(max(abs(scheduleTime - stateTime)) <= tolerance, ...
           'Truncated schedule times do not match the simulated states');

end

function obj = objectiveAtTimeStep(objFunc, simsetup, states, timeStep, ...
                                   model, state, computeStatePartial)

    simsetup.model = model;
    obj = objFunc(simsetup, states, ...
                  'ComputePartials', computeStatePartial, ...
                  'tStep'          , timeStep, ...
                  'state'          , state);

end

%{
  Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
  and SINTEF Digital, Mathematics & Cybernetics.

  This file is part of The Battery Modeling Toolbox BattMo

  BattMo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
%}
