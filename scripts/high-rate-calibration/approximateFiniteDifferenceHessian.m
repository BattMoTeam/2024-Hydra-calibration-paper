function [hessians, report] = approximateFiniteDifferenceHessian( ...
    X, objective, shortnames, varargin)
% Approximate the full objective Hessian by finite differences of gradients.
%
% PerturbationSize may be a positive scalar or vector. Each Hessian is
% computed with respect to the scaled, unit-box parameters. Central
% differences are used where possible. A second-order one-sided stencil is
% used near a bound. The independent Hessian columns and perturbation sizes
% are evaluated in parallel.

    opt = struct('PerturbationSize', 1e-3);
    opt = merge_options(opt, varargin{:});

    X = X(:);
    shortnames = shortnames(:);
    perturbationSizes = opt.PerturbationSize(:)';
    numParams = numel(X);
    numSteps  = numel(perturbationSizes);
    numTasks  = numParams * numSteps;
    columns   = cell(numTasks, 1);
    methods   = cell(numTasks, 1);

    grad0 = evaluateGradient(X, objective, numParams, 'base point');

    parfor idx = 1:numTasks
        [paramidx, stepidx] = ind2sub([numParams, numSteps], idx);

        h = perturbationSizes(stepidx);
        x = X(paramidx);
        direction = zeros(numParams, 1);
        direction(paramidx) = 1;

        doCentral = x - h >= 0 && x + h <= 1;
        doForward = x + 2*h <= 1;
        doBackward = x - 2*h >= 0;

        if doCentral
            method = 'central';
            gradminus = evaluateGradient(X - h*direction, objective, numParams, ...
                                         evaluationLabel(shortnames{paramidx}, h, '-h'));
            gradplus = evaluateGradient(X + h*direction, objective, numParams, ...
                                        evaluationLabel(shortnames{paramidx}, h, '+h'));
            column = (gradplus - gradminus) ./ (2*h);
        elseif doForward
            method = 'forward';
            grad1 = evaluateGradient(X + h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{paramidx}, h, '+h'));
            grad2 = evaluateGradient(X + 2*h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{paramidx}, h, '+2h'));
            column = (-3*grad0 + 4*grad1 - grad2) ./ (2*h);
        elseif doBackward
            method = 'backward';
            grad1 = evaluateGradient(X - h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{paramidx}, h, '-h'));
            grad2 = evaluateGradient(X - 2*h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{paramidx}, h, '-2h'));
            column = (3*grad0 - 4*grad1 + grad2) ./ (2*h);
        else
            error('No feasible Hessian stencil for %s with h=%g.', shortnames{paramidx}, h);
        end

        columns{idx} = column;
        methods{idx} = method;
    end

    allHessians = zeros(numParams, numParams, numSteps);
    allMethods = cell(numParams, numSteps);

    for idx = 1:numTasks
        [paramidx, stepidx] = ind2sub([numParams, numSteps], idx);
        allHessians(:, paramidx, stepidx) = columns{idx};
        allMethods{paramidx, stepidx} = methods{idx};
    end

    hessians = zeros(size(allHessians));
    relAsymmetry = zeros(numSteps, 1);
    relStepChange = NaN(numSteps, 1);

    for stepidx = 1:numSteps
        rawHessian = allHessians(:, :, stepidx);
        hessians(:, :, stepidx) = 0.5 .* (rawHessian + rawHessian');
        relAsymmetry(stepidx) = norm(rawHessian - rawHessian', 'fro') ./ ...
            max(norm(rawHessian, 'fro'), eps);

        if stepidx > 1
            currentHessian = hessians(:, :, stepidx);
            previousHessian = hessians(:, :, stepidx - 1);
            relStepChange(stepidx) = norm(currentHessian - previousHessian, 'fro') ./ ...
                max(norm(currentHessian, 'fro'), eps);
        end
    end

    summary = table(perturbationSizes(:), relAsymmetry, relStepChange, ...
                    'VariableNames', ...
                    {'PerturbationSize', 'RelativeAsymmetry', 'RelativeStepChange'});

    report = struct('allHessians'                , allHessians       , ...
                    'baseGradient'               , baseGradient      , ...
                    'perturbationSizes'          , perturbationSizes , ...
                    'schemes'                    , {allMethods}      , ...
                    'relativeAsymmetry'          , relAsymmetry , ...
                    'relativeStepChange'         , relStepChange, ...
                    'numberOfGradientEvaluations', 1 + 2*numTasks    , ...
                    'summary'                    , summary);

end


function gradient = evaluateGradient(X, objective, numParams, label)

    [objectiveValue, gradient] = objective(X, 'gradientMethod', 'AdjointAD', 'enforceBounds', false);
    gradient = gradient(:);

    assert(isscalar(objectiveValue) && isfinite(objectiveValue), ...
           'Non-finite objective value at %s.', label);
    assert(numel(gradient) == numParams && all(isfinite(gradient)), ...
           'Invalid objective gradient at %s.', label);

end


function label = evaluationLabel(shortname, perturbationSize, offset)

    label = sprintf('%s, h=%g, offset=%s', ...
                    shortname, perturbationSize, offset);

end
