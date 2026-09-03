%% Test finite-difference Hessian approximation on a quadratic objective

clearvars

A = [ 4.0, -1.0,  0.5,  0.0; ...
     -1.0,  3.0,  0.2, -0.4; ...
      0.5,  0.2,  2.0,  0.7; ...
      0.0, -0.4,  0.7,  1.5];
b = [0.3; -0.2; 0.1; 0.4];
objective = @(X, varargin) quadraticObjective(X, A, b, varargin{:});
shortnames = {'p1', 'p2', 'p3', 'p4'};
perturbationSizes = [1e-2, 1e-3, 1e-4];
tolerance = 1e-9;

%% Interior point: central differences

X = [0.2; 0.4; 0.6; 0.8];
[hessians, report] = calculateFDHessian(X, objective, shortnames, perturbationSizes);

for stepIndex = 1:numel(perturbationSizes)
    errorNorm = norm(hessians(:, :, stepIndex) - A, 'fro');
    assert(errorNorm < tolerance, ...
           'Interior Hessian error %g exceeds tolerance for h=%g.', ...
           errorNorm, perturbationSizes(stepIndex));
    assert(all(strcmp(report.schemes(:, stepIndex), 'central')), ...
           'Expected central stencils at the interior point.');
end

%% Bound point: second-order one-sided differences

X = [0; 1; 0.5; 0];
[hessians, report] = calculateFDHessian(X, objective, shortnames, perturbationSizes);

for stepIndex = 1:numel(perturbationSizes)
    errorNorm = norm(hessians(:, :, stepIndex) - A, 'fro');
    assert(errorNorm < tolerance, ...
           'Boundary Hessian error %g exceeds tolerance for h=%g.', ...
           errorNorm, perturbationSizes(stepIndex));
    expectedSchemes = {'forward'; 'backward'; 'central'; 'forward'};
    assert(isequal(report.schemes(:, stepIndex), expectedSchemes), ...
           'Unexpected stencil selection at the bound point.');
end

fprintf('testHessian passed for %d perturbation sizes.\n', ...
        numel(perturbationSizes));


function [value, gradient] = quadraticObjective(X, A, b, varargin) %#ok<INUSD>

    X = X(:);
    value = 0.5 .* X' * A * X + b' * X;
    gradient = A * X + b;

end
