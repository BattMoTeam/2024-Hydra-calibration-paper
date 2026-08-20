function obj = leastSquaresExperimentalVoltage(simsetup, states, expdata, varargin)
%LEASTSQUARESEXPERIMENTALVOLTAGE Voltage mismatch on the simulated domain.
%
% Unlike leastSquaresEI, this objective permits a simulation to terminate
% before all scheduled steps have been completed.  Experimental voltage is
% interpolated only at times for which a nonempty simulated state exists;
% extrapolation is deliberately not enabled.

    opt = struct('ComputePartials', false, ...
                 'tStep'          , []   , ...
                 'state'          , []   , ...
                 'from_states'    , true , ...
                 'maxExtrapolationFraction', 1e-3);
    opt = merge_options(opt, varargin{:});

    states = truncateStatesAtCutoff(states);
    schedule = simsetup.schedule;
    numberOfStates = numel(states);

    assert(numberOfStates > 0, 'The simulation did not produce any states');
    assert(numberOfStates <= numel(schedule.step.val), ...
           'There are more simulated states than scheduled time steps');

    if isempty(opt.tStep)
        timeSteps = (1:numberOfStates)';
    else
        assert(isscalar(opt.tStep) && opt.tStep >= 1 && ...
               opt.tStep <= numberOfStates, ...
               'Requested time step is outside the simulated domain');
        timeSteps = opt.tStep;
    end

    obj = cell(numel(timeSteps), 1);

    for index = 1:numel(timeSteps)
        timeStep = timeSteps(index);
        state = states{timeStep};
        time = state.time;
        dt = schedule.step.val(timeStep);

        lowerOverrun = max(expdata.time(1) - time, 0);
        upperOverrun = max(time - expdata.time(end), 0);
        extrapolation = max(lowerOverrun, upperOverrun);
        roundoffTolerance = 128 * eps(max(abs([time, ...
                                               expdata.time(1), ...
                                               expdata.time(end)])));
        maxExtrapolation = max(opt.maxExtrapolationFraction * dt, ...
                               roundoffTolerance);

        assert(extrapolation <= maxExtrapolation, ...
               ['Simulation time %.17g is outside the experimental domain ', ...
                '[%.17g, %.17g] by %.17g s, which exceeds the allowed ', ...
                'extrapolation of %.17g s'], ...
               time, expdata.time(1), expdata.time(end), extrapolation, ...
               maxExtrapolation);

        if extrapolation > 0
            referenceVoltage = interp1(expdata.time, expdata.U, time, ...
                                       'linear', 'extrap');
        else
            referenceVoltage = interp1(expdata.time, expdata.U, time, 'linear');
        end
        assert(isfinite(referenceVoltage), ...
               'Experimental voltage is unavailable at simulation time %g', time);

        if opt.ComputePartials
            if opt.from_states
                state = simsetup.model.initStateAD(state);
            else
                state = opt.state;
            end
        end

        obj{index} = (state.Control.E - referenceVoltage)^2 * dt;
    end

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
