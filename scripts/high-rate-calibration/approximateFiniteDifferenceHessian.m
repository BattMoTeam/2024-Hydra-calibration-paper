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
    numSteps = numel(perturbationSizes);
    numTasks = numParams * numSteps;
    columns = cell(numTasks, 1);
    taskSchemes = cell(numTasks, 1);

    grad0 = evaluateGradient(X, objective, numParams, 'base point');

    parfor taskIndex = 1:numTasks
        [parameterIndex, stepIndex] = ind2sub([numParams, numSteps], taskIndex);

        h = perturbationSizes(stepIndex);
        x = X(parameterIndex);
        direction = zeros(numParams, 1);
        direction(parameterIndex) = 1;

        okCentral = x - h >= 0 && x + h <= 1;
        okForward = x + 2*h <= 1;
        okBackward = x - 2*h >= 0;

        if okCentral
            scheme = 'central';
            gradminus = evaluateGradient(X - h*direction, objective, numParams, ...
                                         evaluationLabel(shortnames{parameterIndex}, h, '-h'));
            gradplus = evaluateGradient(X + h*direction, objective, numParams, ...
                                        evaluationLabel(shortnames{parameterIndex}, h, '+h'));
            column = (gradplus - gradminus) ./ (2*h);
        elseif okForward
            scheme = 'forward';
            grad1 = evaluateGradient(X + h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{parameterIndex}, h, '+h'));
            grad2 = evaluateGradient(X + 2*h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{parameterIndex}, h, '+2h'));
            column = (-3*grad0 + 4*grad1 - grad2) ./ (2*h);
        elseif okBackward
            scheme = 'backward';
            grad1 = evaluateGradient(X - h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{parameterIndex}, h, '-h'));
            grad2 = evaluateGradient(X - 2*h*direction, objective, numParams, ...
                                     evaluationLabel(shortnames{parameterIndex}, h, '-2h'));
            column = (3*grad0 - 4*grad1 + grad2) ./ (2*h);
        else
            error('No feasible Hessian stencil for %s with h=%g.', ...
                  shortnames{parameterIndex}, h);
        end

        columns{taskIndex} = column;
        taskSchemes{taskIndex} = scheme;
    end

    rawHessians = zeros(numParams, numParams, numSteps);
    schemes = cell(numParams, numSteps);

    for taskIndex = 1:numTasks
        [parameterIndex, stepIndex] = ind2sub([numParams, numSteps], taskIndex);
        rawHessians(:, parameterIndex, stepIndex) = columns{taskIndex};
        schemes{parameterIndex, stepIndex} = taskSchemes{taskIndex};
    end

    hessians = zeros(size(rawHessians));
    relAsymmetry = zeros(numSteps, 1);
    relStepChange = NaN(numSteps, 1);

    for stepIndex = 1:numSteps
        rawHessian = rawHessians(:, :, stepIndex);
        hessians(:, :, stepIndex) = 0.5 .* (rawHessian + rawHessian');
        relAsymmetry(stepIndex) = norm(rawHessian - rawHessian', 'fro') ./ ...
            max(norm(rawHessian, 'fro'), eps);

        if stepIndex > 1
            currentHessian = hessians(:, :, stepIndex);
            previousHessian = hessians(:, :, stepIndex - 1);
            relStepChange(stepIndex) = norm(currentHessian - previousHessian, 'fro') ./ ...
                max(norm(currentHessian, 'fro'), eps);
        end
    end

    summary = table(perturbationSizes(:), relAsymmetry, relStepChange, ...
                    'VariableNames', ...
                    {'PerturbationSize', 'RelativeAsymmetry', 'RelativeStepChange'});

    report = struct('rawHessians'                , rawHessians       , ...
                    'baseGradient'               , baseGradient      , ...
                    'perturbationSizes'          , perturbationSizes , ...
                    'schemes'                    , {schemes}         , ...
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
