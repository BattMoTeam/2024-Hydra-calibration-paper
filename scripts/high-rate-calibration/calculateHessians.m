function [Hscaled, HfdScaled] = calculateHessians( ...
        invHscaled, Xopt, objective, shortnames, debug, ...
        hessianSteps, hessianfdpertsize)
% Calculate and compare the BFGS and finite-difference Hessians.

    issym = @(x) max(max(abs(x-x'))) < 1e-10;
    symmetrize = @(x) 0.5.*(x + x');

    assert(issym(invHscaled), 'invHscaled is not symmetric');
    invHscaled = symmetrize(invHscaled);
    Hscaled = invHscaled \ eye(numel(Xopt));
    assert(issym(Hscaled), 'Hscaled is not symmetric');
    Hscaled = symmetrize(Hscaled);

    [eigenvecsBFGS, eigenvalsBFGS] = eig(Hscaled, 'vector');
    disp(shortnames);

    for k = 1:numel(Xopt)
        fprintf('Eigenvalue  %d: %g\n', k, eigenvalsBFGS(k));
        fprintf('Eigenvector %d: ', k);
        fprintf('%g ', eigenvecsBFGS(:, k));
        fprintf('\n');
    end

    % Fix signs: make the largest absolute value in each eigenvector positive
    for k = 1:numel(Xopt)
        [~, maxIdx] = max(abs(eigenvecsBFGS(:, k)));
        if eigenvecsBFGS(maxIdx, k) < 0
            eigenvecsBFGS(:, k) = -eigenvecsBFGS(:, k);
        end
    end
    disp('After fixing signs:');
    for k = 1:numel(Xopt)
        fprintf('Eigenvalue  %d: %g\n', k, eigenvalsBFGS(k));
        fprintf('Eigenvector %d: ', k);
        fprintf('%g ', eigenvecsBFGS(:, k));
        fprintf('\n');
    end

    [U, singularvals, V] = svd(Hscaled, 'econ'); %#ok<ASGLU>
    singularvals = diag(singularvals);

    reltol = 1e-10;
    ranktol = reltol * singularvals(1);
    numericalRank = nnz(singularvals > ranktol);
    condno = singularvals(1) / singularvals(end);

    fprintf('Numerical rank: %d/%d\n', numericalRank, numel(Xopt));
    fprintf('SVD condition number: %.3e\n', condno);
    fprintf('Negative eigenvals: %d\n', nnz(eigenvalsBFGS < -ranktol));

    if debug
        [HfdScaled, HfdReport] = approximateFiniteDifferenceHessian( ...
            Xopt, objective, shortnames, 'PerturbationSize', hessianSteps);

        numSteps = numel(hessianSteps);
        relErr = zeros(numSteps, 1);
        maxAbsErr = zeros(numSteps, 1);
        fdeigenvals = zeros(numel(Xopt), numSteps);

        for stepno = 1:numSteps
            Hfd = HfdScaled(:, :, stepno);
            Hdiff = Hscaled - Hfd;

            relErr(stepno) = norm(Hdiff, 'fro') ./ max(norm(Hfd, 'fro'), eps);
            maxAbsErr(stepno) = max(abs(Hdiff), [], 'all');
            fdeigenvals(:, stepno) = sort(eig(Hfd));
        end

        Hcomparison = HfdReport.summary;
        Hcomparison.relErrToBFGS = relErr;
        Hcomparison.maxAbsErr = maxAbsErr;
        disp('Finite-difference Hessian comparison:');
        disp(Hcomparison);
        fprintf('Gradient evaluations for finite-difference Hessians: %d\n', ...
                HfdReport.numberOfGradientEvaluations);

        stencilNames = matlab.lang.makeUniqueStrings( ...
            matlab.lang.makeValidName(compose('h_%g', hessianSteps)));
        stencilTable = cell2table(HfdReport.schemes, ...
                                  'VariableNames', cellstr(stencilNames), ...
                                  'RowNames', shortnames);
        disp('Finite-difference stencil by parameter:');
        disp(stencilTable);

        eigenvalueComparison = table((1:numel(Xopt))', eigenvalsBFGS, ...
                                     'VariableNames', {'Mode', 'BFGS'});
        for stepno = 1:numSteps
            varname = matlab.lang.makeValidName( ...
                sprintf('FD_h_%g', hessianSteps(stepno)));
            eigenvalueComparison.(varname) = fdeigenvals(:, stepno);
        end
        disp('Hessian eigenvalue comparison:');
        disp(eigenvalueComparison);

        boundtol = 1e-8;
        freeParams = Xopt > boundtol & Xopt < 1 - boundtol;
        activeParams = ~freeParams;

        if any(activeParams)
            activeParameterTable = table( ...
                shortnames(activeParams), Xopt(activeParams), ...
                'VariableNames', {'Parameter', 'ScaledValue'});
            disp('Parameters active at a unit-box bound:');
            disp(activeParameterTable);
        end

        fprintf('Free parameters for reduced-Hessian analysis: %d/%d\n', ...
                nnz(freeParams), numel(Xopt));
        if any(freeParams)
            for stepno = 1:numSteps
                Hreduced = HfdScaled(freeParams, freeParams, stepno);
                reducedEigenvalues = eig(Hreduced);
                fprintf(['FD reduced Hessian h=%g: min eigenvalue=%g, ', ...
                         'max eigenvalue=%g\n'], ...
                        hessianSteps(stepno), ...
                        min(reducedEigenvalues), max(reducedEigenvalues));
            end
        end
    end

    [HfdScaled, HfdReport] = approximateFiniteDifferenceHessian( ... %#ok<ASGLU>
        Xopt, objective, shortnames, ...
        'PerturbationSize', hessianfdpertsize);

end
