function l2 = l2error(x1, y1, x2, y2, varargin)

    assert(numel(x1) == numel(y1), 'x1 and y1 must have the same number of elements');
    assert(numel(x2) == numel(y2), 'x2 and y2 must have the same number of elements');

    opt = struct('weight'  , true , ...
                 'extrap'  , false, ...
                 'truncate', false);
    opt = merge_options(opt, varargin{:});

    % Truncate (x2, y2) to the range of x1
    if opt.truncate && x2(end) > x1(end)
        idx = find(x2 > x1(end), 1, 'first');
        x2 = x2(1:idx);
        y2 = y2(1:idx);
    end

    % Interpolate to x1
    if opt.extrap
        y2i = interp1(x2, y2, x1, 'linear', 'extrap');
    else
        y2i = interp1(x2, y2, x1);
    end

    assert(all(isfinite(y2i)), 'Interpolated values are not finite. Check if x1 is within the range of x2 (both start and end)');

    % Compute l2 error
    l2sq = trapz(x1, (y1 - y2i).^2);

    if opt.weight
        l2sq = l2sq / (x1(end) - x1(1));
    end

    l2 = sqrt(l2sq);

end
