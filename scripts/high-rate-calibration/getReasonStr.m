function [reasonStr, tbl] = getReasonStr(history, colname)

    if nargin < 2
        colname = '';
    end

    if numel(history.val) == 1

        warning('Only one iteration in callOptimizer');
        valstr = sprintf('Value is %g\n', history.val(end));
        valdiffstr = 'No value diff\n';
        pgstr = sprintf('Gradient value is %g\n', history.pg(end));

    elseif numel(history.val) == 2

        warning('Only two iterations in callOptimizer');
        diffop = @(v, offset) abs(v(end)-v(end-1));
        valstr = sprintf('Value is %g\n', history.val(end));
        valdiffstr = sprintf('Value diff is %g\n', diffop(history.val, 0));
        pgstr = sprintf('Gradient diff is %g\n', diffop(history.pg, 0));

    else

        diffop = @(v, offset) abs(v(end-offset)-v(end-offset-1));
        valstr = sprintf('Value is %g\n', history.val(end));
        valdiffstr = sprintf('Value diffs (prev last %g) %g\n', diffop(history.val, 1), diffop(history.val, 0));
        pgstr = sprintf('Gradient diffs (prev last %g) %g\n', diffop(history.pg, 1), diffop(history.pg, 0));

    end

    reasonStr = [sprintf('Reason for termination:\n'), ...
                 valstr, ...
                 valdiffstr, ...
                 pgstr, ...
                 sprintf('number of iterations %g\n', numel(history.val))];

    rows = {'Obj value (end)'; 'Obj value diff (end-1:end)'; 'Pg (end)'; 'Pg diff (end-1:end)'; 'Num iterations'};
    values = [history.val(end), ...
              abs(history.val(end)-history.val(end-1)), ...
              history.pg(end), ...
              abs(history.pg(end)-history.pg(end-1)), ...
              numel(history.val)]';

    tbl = table(values, 'RowNames', rows);

    if ~isempty(colname)
        tbl.Properties.VariableNames = {colname};
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

  BattMo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
