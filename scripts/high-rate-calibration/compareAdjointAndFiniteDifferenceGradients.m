function [tbl, comparison] = compareAdjointAndFiniteDifferenceGradients( ...
        X, objective, shortnames, varargin)
% Compare adjoint and finite-difference gradients.
%
% PerturbationSize may be a positive scalar or vector and defaults to
% 1e-5. A vector requests one finite-difference comparison per step size.
% adjointGradient may be supplied to avoid recomputing the adjoint.
% DifferenceScheme may be 'central' (default) or 'forward'. Central
% differences fall back to a feasible one-sided difference at a bound.

    opt = struct( ...
        'PerturbationSize', 1e-5, ...
        'adjointGradient', [], ...
        'DifferenceScheme', 'central');
    opt = merge_options(opt, varargin{:});

    differenceScheme = validatestring( ...
        opt.DifferenceScheme, {'central', 'forward'}, ...
        mfilename, 'DifferenceScheme');

    perturbationSizes = opt.PerturbationSize;
    validateattributes( ...
        perturbationSizes, {'numeric'}, ...
        {'real', 'finite', 'positive', 'nonempty', 'vector', '<=', 1}, ...
        mfilename, 'PerturbationSize');
    perturbationSizes = perturbationSizes(:)';

    if isempty(opt.adjointGradient)
        [~, gad] = objective(X, 'gradientMethod', 'AdjointAD');
    else
        gad = opt.adjointGradient;
    end

    gad = gad(:);
    shortnames = shortnames(:);

    nPert = numel(perturbationSizes);
    nParam = numel(X);
    finiteDifferenceGradients = cell(nPert, 1);
    schemes = cell(nParam, nPert);

    assert(numel(shortnames) == nParam, ...
           'The shortname and parameter lists must have the same length.');
    assert(numel(gad) == nParam, ...
           'The adjoint gradient and parameter lists must have the same length.');
    assert(all(X(:) >= 0 & X(:) <= 1), ...
           'Scaled parameters must lie in the interval [0, 1].');

    for k = 1:nPert
        h = perturbationSizes(k);

        firstSteps = cell(nParam, 1);
        secondSteps = cell(nParam, 1);

        for ip = 1:nParam
            canStepForward = X(ip) + h <= 1;
            canStepBackward = X(ip) - h >= 0;
            assert(canStepForward || canStepBackward, ...
                   ['Perturbation size %g has no feasible direction for ', ...
                    'scaled parameter %s.'], h, shortnames{ip});

            if strcmp(differenceScheme, 'central') && ...
                    canStepForward && canStepBackward
                firstSteps{ip} = h;
                secondSteps{ip} = -h;
                schemes{ip, k} = 'central';
            elseif canStepForward
                firstSteps{ip} = h;
                secondSteps{ip} = h;
                schemes{ip, k} = 'forward';
            else
                firstSteps{ip} = -h;
                secondSteps{ip} = -h;
                schemes{ip, k} = 'backward';
            end
        end

        [~, firstGradient] = objective( ...
            X, ...
            'gradientMethod', 'PerturbationADNUM', ...
            'PerturbationSize', firstSteps);

        if strcmp(differenceScheme, 'central')
            [~, secondGradient] = objective( ...
                X, ...
                'gradientMethod', 'PerturbationADNUM', ...
                'PerturbationSize', secondSteps);
            finiteDifferenceGradients{k} = ...
                0.5 .* (firstGradient + secondGradient);
        else
            finiteDifferenceGradients{k} = firstGradient;
        end

        finiteDifferenceGradients{k} = finiteDifferenceGradients{k}(:);
    end

    tbl = table(shortnames, gad, ...
                'VariableNames', {'Shortname', 'Adjoint'});

    finiteDifferenceGradient = zeros(nParam, nPert);
    absoluteError = zeros(nParam, nPert);
    relativeError = zeros(nParam, nPert);

    for k = 1:nPert
        h = perturbationSizes(k);
        gnum = finiteDifferenceGradients{k};
        absErr = abs(gad - gnum);
        denominator = max(abs(gad), abs(gnum));
        relErr = absErr ./ denominator;

        % Relative error is meaningless when both gradients are tiny.
        relErr(denominator < 1e-10) = NaN;

        finiteDifferenceGradient(:, k) = gnum;
        absoluteError(:, k) = absErr;
        relativeError(:, k) = relErr;

        suffix = matlab.lang.makeValidName(sprintf('h_%g', h));
        tbl.(['FD_' suffix]) = gnum;
        tbl.(['AbsErr_' suffix]) = absErr;
        tbl.(['RelErr_' suffix]) = relErr;
    end

    relativeTolerance = 1e-2;
    comparison = struct( ...
        'adjointGradient', gad, ...
        'finiteDifferenceGradient', finiteDifferenceGradient, ...
        'absoluteError', absoluteError, ...
        'relativeError', relativeError, ...
        'passed', relativeError <= relativeTolerance, ...
        'perturbationSize', perturbationSizes, ...
        'scheme', {schemes}, ...
        'requestedScheme', differenceScheme, ...
        'relativeTolerance', relativeTolerance);

    disp(tbl);

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
