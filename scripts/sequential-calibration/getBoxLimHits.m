function tblOpt = getBoxLimHits(simSetup, params, PS, Xopt)

    lims = params{1}.boxLims;
    names = PS.shortnames;

    % Get optimal params
    setupOpt = updateSetupFromScaledParameters(simSetup, params, Xopt);
    pOpt = applyFunction(@(p) p.getParameter(setupOpt), params);
    pOpt = vertcat(pOpt{:});

    % Rel error
    p0 = applyFunction(@(p) p.getParameter(simSetup), params);
    p0 = vertcat(p0{:});
    relchange = abs(p0 - pOpt) ./ abs(p0);

    % Check boxLim hits
    lower = lims(:, 1);
    upper = lims(:, 2);
    tol = 1e-3;
    relErrLower = abs(pOpt - lower) ./ abs(pOpt);
    relErrUpper = abs(pOpt - upper) ./ abs(pOpt);
    assert(all(isfinite(relErrLower) & isfinite(relErrUpper)), 'Relative errors must be finite. Check if pOpt is zero or if lims are zero.');
    isHitLower = relErrLower < tol;
    isHitUpper = relErrUpper < tol;
    isHit = isHitLower | isHitUpper;
    boxLimHit = char(isHit * '*');

    if isscalar(p0)
        % One-row table can't have character vector, so convert to string
        boxLimHit = string(boxLimHit);
    end

    % Print as a table
    tblOpt = table(p0, pOpt, relchange, boxLimHit, 'rownames', names);

    if any(isHit)
        warning('boxLim is hit. Check boxLims.');
    end

end
