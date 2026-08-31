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
    perturbationSizes = opt.PerturbationSize;

    validateattributes( ...
        perturbationSizes, {'numeric'}, ...
        {'real', 'finite', 'positive', 'nonempty', 'vector', '<=', 0.25}, ...
        mfilename, 'PerturbationSize');
    perturbationSizes = perturbationSizes(:)';

    numberOfParameters = numel(X);
    numberOfStepSizes = numel(perturbationSizes);

    assert(numel(shortnames) == numberOfParameters, ...
           'The shortname and parameter lists must have the same length.');
    assert(all(isfinite(X)) && all(X >= 0 & X <= 1), ...
           'Scaled parameters must be finite and lie in the interval [0, 1].');

    baseGradient = evaluateGradient(X, objective, numberOfParameters, 'base point');

    numberOfTasks = numberOfParameters * numberOfStepSizes;
    columns = cell(numberOfTasks, 1);
    taskSchemes = cell(numberOfTasks, 1);

    parfor taskIndex = 1:numberOfTasks
        [parameterIndex, stepIndex] = ind2sub( ...
            [numberOfParameters, numberOfStepSizes], taskIndex);

        h = perturbationSizes(stepIndex);
        x = X(parameterIndex);
        direction = zeros(numberOfParameters, 1);
        direction(parameterIndex) = 1;

        canStepCentral = x - h >= 0 && x + h <= 1;
        canStepForward = x + 2*h <= 1;
        canStepBackward = x - 2*h >= 0;

        if canStepCentral
            scheme = 'central';
            gradientMinus = evaluateGradient( ...
                X - h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '-h'));
            gradientPlus = evaluateGradient( ...
                X + h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '+h'));
            column = (gradientPlus - gradientMinus) ./ (2*h);
        elseif canStepForward
            scheme = 'forward';
            gradientOne = evaluateGradient( ...
                X + h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '+h'));
            gradientTwo = evaluateGradient( ...
                X + 2*h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '+2h'));
            column = (-3*baseGradient + 4*gradientOne - gradientTwo) ./ (2*h);
        elseif canStepBackward
            scheme = 'backward';
            gradientOne = evaluateGradient( ...
                X - h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '-h'));
            gradientTwo = evaluateGradient( ...
                X - 2*h*direction, objective, numberOfParameters, ...
                evaluationLabel(shortnames{parameterIndex}, h, '-2h'));
            column = (3*baseGradient - 4*gradientOne + gradientTwo) ./ (2*h);
        else
            error('No feasible Hessian stencil for %s with h=%g.', ...
                  shortnames{parameterIndex}, h);
        end

        columns{taskIndex} = column;
        taskSchemes{taskIndex} = scheme;
    end

    rawHessians = zeros(numberOfParameters, numberOfParameters, numberOfStepSizes);
    schemes = cell(numberOfParameters, numberOfStepSizes);

    for taskIndex = 1:numberOfTasks
        [parameterIndex, stepIndex] = ind2sub( ...
            [numberOfParameters, numberOfStepSizes], taskIndex);
        rawHessians(:, parameterIndex, stepIndex) = columns{taskIndex};
        schemes{parameterIndex, stepIndex} = taskSchemes{taskIndex};
    end

    hessians = zeros(size(rawHessians));
    relativeAsymmetry = zeros(numberOfStepSizes, 1);
    relativeStepChange = NaN(numberOfStepSizes, 1);

    for stepIndex = 1:numberOfStepSizes
        rawHessian = rawHessians(:, :, stepIndex);
        hessians(:, :, stepIndex) = 0.5 .* (rawHessian + rawHessian');
        relativeAsymmetry(stepIndex) = ...
            norm(rawHessian - rawHessian', 'fro') ./ ...
            max(norm(rawHessian, 'fro'), eps);

        if stepIndex > 1
            currentHessian = hessians(:, :, stepIndex);
            previousHessian = hessians(:, :, stepIndex - 1);
            relativeStepChange(stepIndex) = ...
                norm(currentHessian - previousHessian, 'fro') ./ ...
                max(norm(currentHessian, 'fro'), eps);
        end
    end

    summary = table( ...
        perturbationSizes(:), relativeAsymmetry, relativeStepChange, ...
        'VariableNames', ...
        {'PerturbationSize', 'RelativeAsymmetry', 'RelativeStepChange'});

    report = struct( ...
        'rawHessians'               , rawHessians, ...
        'baseGradient'              , baseGradient, ...
        'perturbationSizes'         , perturbationSizes, ...
        'schemes'                   , {schemes}, ...
        'relativeAsymmetry'         , relativeAsymmetry, ...
        'relativeStepChange'        , relativeStepChange, ...
        'numberOfGradientEvaluations', 1 + 2*numberOfTasks, ...
        'summary'                   , summary);

end


function gradient = evaluateGradient(X, objective, numberOfParameters, label)

    [objectiveValue, gradient] = objective( ...
        X, 'gradientMethod', 'AdjointAD', 'enforceBounds', false);
    gradient = gradient(:);

    assert(isscalar(objectiveValue) && isfinite(objectiveValue), ...
           'Non-finite objective value at %s.', label);
    assert(numel(gradient) == numberOfParameters && all(isfinite(gradient)), ...
           'Invalid objective gradient at %s.', label);

end


function label = evaluationLabel(shortname, perturbationSize, offset)

    label = sprintf('%s, h=%g, offset=%s', ...
                    shortname, perturbationSize, offset);

end
