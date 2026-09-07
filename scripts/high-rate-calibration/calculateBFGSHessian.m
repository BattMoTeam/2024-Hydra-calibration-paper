function Hscaled = calculateBFGSHessian(invHscaled, shortnames)
% Calculate the BFGS Hessian from the inverse Hessian stored in the optimization history.

    relnorm = @(x) norm(x-x', 'fro') ./ max(norm(x, 'fro'), eps);
    issym = @(x) relnorm(x) < 1e-10;
    symmetrize = @(x) 0.5 .* (x + x');

    assert(issym(invHscaled), 'invHscaled is not symmetric %g', relnorm(invHscaled));
    invHscaled = symmetrize(invHscaled);
    Hscaled = invHscaled \ eye(size(invHscaled));
    assert(issym(Hscaled), 'Hscaled is not symmetric: %g', relnorm(Hscaled));
    Hscaled = symmetrize(Hscaled);

    [eigenvecs, eigenvals] = eig(Hscaled, 'vector');
    disp(shortnames);

    % disp('BFGS H:');
    % for k = 1:numel(eigenvals)
    %     fprintf('Eigenvalue  %d: %g\n', k, eigenvals(k));
    %     fprintf('Eigenvector %d: ', k);
    %     fprintf('%g ', eigenvecs(:, k));
    %     fprintf('\n');
    % end

    % Fix signs: make the largest absolute value in each eigenvector positive
    for k = 1:numel(eigenvals)
        [~, maxIdx] = max(abs(eigenvecs(:, k)));
        if eigenvecs(maxIdx, k) < 0
            eigenvecs(:, k) = -eigenvecs(:, k);
        end
    end

    disp('BFGS H After fixing signs:');
    for k = 1:numel(eigenvals)
        fprintf('Eigenvalue  %d: %g\n', k, eigenvals(k));
        fprintf('Eigenvector %d: ', k);
        fprintf('%g ', eigenvecs(:, k));
        fprintf('\n');
    end

    [U, singularvals, V] = svd(Hscaled, 'econ'); %#ok<ASGLU>
    singularvals = diag(singularvals);

    reltol = 1e-10;
    ranktol = reltol * singularvals(1);
    numericalRank = nnz(singularvals > ranktol);
    condno = singularvals(1) / singularvals(end);

    fprintf('BFGS H Numerical rank: %d/%d\n', numericalRank, size(Hscaled, 1));
    fprintf('BFGS H SVD condition number: %.3e\n', condno);

end
