function tbl = compareAdjointAndFiniteDifferenceGradients(X, objective, shortnames)

    perturbationSizes = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];

    [~, gad] = objective(X, 'gradientMethod', 'AdjointAD');

    gad = gad(:);
    shortnames = shortnames(:);

    nPert = numel(perturbationSizes);
    nParam = numel(X);
    finiteDifferenceGradients = cell(nPert, 1);

    assert(numel(shortnames) == nParam, ...
           'The shortname and parameter lists must have the same length.');

    for k = 1:nPert
        h = perturbationSizes(k);

        % Perturb in a feasible direction for each scaled parameter.
        parameterSteps = repmat({h}, nParam, 1);

        for ip = 1:nParam
            if X(ip) + h > 1
                parameterSteps{ip} = -h;
            end
        end

        [~, finiteDifferenceGradients{k}] = objective( ...
            X, ...
            'gradientMethod', 'PerturbationADNUM', ...
            'PerturbationSize', parameterSteps);

        finiteDifferenceGradients{k} = finiteDifferenceGradients{k}(:);
    end

    tbl = table(shortnames, gad, ...
                'VariableNames', {'Shortname', 'Adjoint'});

    for k = 1:nPert
        h = perturbationSizes(k);
        gnum = finiteDifferenceGradients{k};
        absErr = abs(gad - gnum);
        denominator = max(abs(gad), abs(gnum));
        relErr = absErr ./ denominator;

        % Relative error is meaningless when both gradients are tiny.
        relErr(denominator < 1e-10) = NaN;

        suffix = matlab.lang.makeValidName(sprintf('h_%g', h));
        tbl.(['FD_' suffix]) = gnum;
        tbl.(['AbsErr_' suffix]) = absErr;
        tbl.(['RelErr_' suffix]) = relErr;
    end

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
