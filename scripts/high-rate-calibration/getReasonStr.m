function [reasonStr, tbl] = getReasonStr(history, varargin)

    opt = struct('gradTol'     , [], ...
                 'objChangeTol', [], ...
                 'maxit'       , []);
    opt = merge_options(opt, varargin{:});

    % Pass all or none
    criteriaPassed = [~isempty(opt.gradTol), ...
                      ~isempty(opt.objChangeTol), ...
                      ~isempty(opt.maxit)];
    assert(all(criteriaPassed) || ~any(criteriaPassed), ...
           'Pass gradTol, objChangeTol, and maxit together.');
    markCriteria = all(criteriaPassed);

    if numel(history.val) == 1

        warning('Only one iteration in callOptimizer');
        valstr = sprintf('Value is %g\n', history.val(end));
        valdiffstr = 'No value diff\n';
        pgstr = sprintf('Pg value is %g\n', history.pg(end));
        pgdiffstr = 'No pg diff\n';

    elseif numel(history.val) == 2

        warning('Only two iterations in callOptimizer');
        diffop = @(v, offset) abs(v(end)-v(end-1));
        valstr = sprintf('Value is %g\n', history.val(end));
        valdiffstr = sprintf('Value diff is %g\n', diffop(history.val, 0));
        pgstr = sprintf('Pg value is %g\n', history.pg(end));
        pgdiffstr = sprintf('Pg diff is %g\n', diffop(history.pg, 0));

    else

        diffop = @(v, offset) abs(v(end-offset)-v(end-offset-1));
        valstr = sprintf('Values (prev last %g) %g\n', history.val(end-1), history.val(end));
        valdiffstr = sprintf('Value diffs (prev last %g) %g\n', diffop(history.val, 1), diffop(history.val, 0));
        pgstr = sprintf('Pg values (prev last %g) %g\n', history.pg(end-1), history.pg(end));
        pgdiffstr = sprintf('Pg diffs (prev last %g) %g\n', diffop(history.pg, 1), diffop(history.pg, 0));

    end

    numIterations = numel(history.val);
    % The history also contains the initial point at BFGS iteration zero.
    optimizerIterations = numIterations - 1;
    itstr = sprintf('number of iterations %g\n', numIterations);

    % Put a > in front of the criteria (may be several)
    if markCriteria
        if numel(history.val) >= 2 && ...
                abs(history.val(end) - history.val(end-1)) < opt.objChangeTol
            valdiffstr = ['> ', valdiffstr];
        end
        if history.pg(end) < opt.gradTol
            pgstr = ['> ', pgstr];
        end
        if optimizerIterations >= opt.maxit
            itstr = ['> ', itstr];
        end
    end

    reasonStr = [sprintf('Reason for termination:\n'), ...
                 valstr, ...
                 valdiffstr, ...
                 pgstr, ...
                 pgdiffstr, ...
                 itstr];

    rows = {'Obj value (end)'; 'Obj value diff (end-1:end)'; 'Pg (end)'; 'Pg diff (end-1:end)'; 'Num iterations'};
    values = [history.val(end), ...
              abs(history.val(end)-history.val(end-1)), ...
              history.pg(end), ...
              abs(history.pg(end)-history.pg(end-1)), ...
              numIterations]';

    tbl = table(values, 'RowNames', rows);

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
