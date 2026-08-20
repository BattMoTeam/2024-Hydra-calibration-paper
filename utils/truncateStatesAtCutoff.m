function states = truncateStatesAtCutoff(states, cutoffVoltage)
%TRUNCATESTATESATCUTOFF Remove empty states and data beyond voltage cutoff.
%
% Packed simulations can return a cell array with unused entries after an
% early stop.  Remove those entries before consumers use cellfun or index
% the state sequence.  If a cutoff is supplied, retain only the states
% strictly above that cutoff.  The first state at or below cutoff is used
% by the forward solver to detect termination, but is not part of the
% comparison domain or adjoint schedule.

    if nargin < 2
        cutoffVoltage = [];
    end

    assert(iscell(states), 'states must be a cell array');
    states = states(~cellfun(@isempty, states));

    if isempty(states) || isempty(cutoffVoltage)
        return
    end

    voltage = cellfun(@(state) state.Control.E, states);
    cutoffIndex = find(voltage <= cutoffVoltage, 1, 'first');

    if ~isempty(cutoffIndex)
        states = states(1:(cutoffIndex - 1));
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
