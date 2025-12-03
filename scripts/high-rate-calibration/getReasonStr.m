function reasonStr = getReasonStr(history)

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

end
