function reason = getBFGSstopReason(history, bfgsopts0)

    hasbfgsopts = nargin >= 2;

    if hasbfgsopts
        if isstruct(bfgsopts0)
            bfgsopts = bfgsopts0;
        elseif iscell(bfgsopts0)
            names  = bfgsopts0(1:2:end);
            values = bfgsopts0(2:2:end);
            bfgsopts = struct();
            for k = 1:numel(names)
                bfgsopts.(names{k}) = values{k};
            end
        end
    end

    diffop = @(v, offset) abs(v(end)-v(end-1));

    if isscalar(history.val)

        warning('Only one iteration in history');
        valstr = sprintf('Obj val is %g\n', history.val(end));
        pgstr = sprintf('Pg value is %g\n', history.pg(end));
        alphastr = sprintf('Step size is %g\n', history.alpha(end));

    elseif numel(history.val) == 2

        warning('Only two iterations in history');
        valstr = sprintf('Obj vals are %g and %g (diff %g)\n', history.val(end-1), history.val(end), diffop(history.val, 0));
        pgstr = sprintf('Pg vals are %g and %g (diff %g)\n', history.pg(end-1), history.pg(end), diffop(history.pg, 0));
        alphastr = sprintf('Step sizes are %g and %g\n', history.alpha(end-1), history.alpha(end));

        if hasbfgsopts
            if diffop(history.val, 0) <= bfgsopts.objChangeTol
                valstr = ['-> ', valstr];
            end
        end

    else

        valstr = sprintf('Obj vals are [..., %g, %g, %g] (diffs %g, %g)\n', history.val(end-2), history.val(end-1), history.val(end), diffop(history.val, 1), diffop(history.val, 0));
        pgstr = sprintf('Pg vals are [..., %g, %g, %g] (diffs %g, %g)\n', history.pg(end-2), history.pg(end-1), history.pg(end), diffop(history.pg, 1), diffop(history.pg, 0));
        alphastr = sprintf('Step sizes are [..., %g, %g, %g]\n', history.alpha(end-2), history.alpha(end-1), history.alpha(end));

        if hasbfgsopts
            if diffop(history.val, 0) <= bfgsopts.objChangeTol
                valstr = ['-> ', valstr];
            end
        end

    end

    itstr = sprintf('number of iterations %g (maxit=%g)\n', numel(history.val), bfgsopts.maxit);

    if hasbfgsopts
        if history.pg(end) <= bfgsopts.gradTol
            pgstr = ['-> ', pgstr];
        end

        if numel(history.val) >= bfgsopts.maxit
            itstr = ['-> ', itstr];
        end
    end

    reason = [sprintf('Reason for BFGS stopping:\n'), ...
              valstr, ...
              pgstr, ...
              alphastr, ...
              itstr];

end
