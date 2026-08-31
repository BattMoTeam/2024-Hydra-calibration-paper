%% Script to calibrate parameters using high-rate data

clear all
close all

diaryname = sprintf('_diary-%s-%s.txt', mfilename, datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
diary(diaryname);

am    = 'ActiveMaterial';
itf   = 'Interface';
pe    = 'PositiveElectrode';
ne    = 'NegativeElectrode';
co    = 'Coating';
sd    = 'SolidDiffusion';
ctrl  = 'Control';
elyte = 'Electrolyte';
sep   = 'Separator';

% mrstDebug(0);

doplot = true;
debug = true;
hessian = true;
dosave = true;
gradientStepSizes = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
hessianStepSizes = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9];

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

%% Fetch experimental data

datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Highest DRate is last
k = numel(dataraw.time);
expdata = struct('time', dataraw.time{k} * hour, ...
                 'U'   , dataraw.voltage{k}    , ...
                 'I'   , abs(mean(dataraw.current{k})));

%% Initial guess using equilibrium calibration data

filename     = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
jsonstructEC = parseBattmoJson(filename);

%shortnames = {'ne_vsa', 'pe_vsa', 'ne_bg', 'pe_bg', 'ne_D', 'pe_D', 'elyte_bgfactorKappa', 'elyte_bgfactorD'};
% Use 'elyte_bgfactor' instead of the two specific factors to calibrate one shared factor.
% shortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', 'elyte_bgfactorKappa', 'elyte_bgfactorD'};
shortnames = {'pe_vsa', 'ne_D', 'pe_D', 'elyte_bgfactor'};
% shortnames = {'pe_vsa', 'ne_D', 'pe_D', 'elyte_bgfactorKappa', 'elyte_bgfactorD'};
% shortnames = {'pe_vsa', 'ne_D', 'pe_D', 'elyte_bgfactorD'};
disp('shortnames:');
printer(shortnames);
useRegionBruggemanCoefficients = any(contains(shortnames, 'elyte_bg'));

numTimesteps = 100; % 400

input0 = struct('I'                             , expdata.I, ...
                'totalTime'                     , expdata.time(end)             , ...
                'numTimesteps'                  , numTimesteps                  , ...
                'lowRateParams'                 , jsonstructEC                  , ...
                'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                'include_current_collectors'    , true);
output0 = runHydra(input0, 'clearSimulation', false);

% Avoid adaptive time stepping to ensure adjoint consistency
output0.nls.timeStepSelector = SimpleTimeStepSelector();
output0.nls.maxTimestepCuts = 0;

if debug
    % Check how exp and initial guess compare
    figure; hold on; grid on;
    plot(expdata.time/hour, expdata.U, 'k--');
    plot(getTime(output0.states)/hour, getE(output0.states));
    xlabel('time / h')
    ylabel('potential / V')
    title('initial guess')
    drawnow
end

%% Setup optimization

% Evaluate experimental data at simulation times (allow for
% extrapolation since expdata.time(end) is very close to
% output.states{end}.time)
simtimes = getTime(output0.states);
assert(expdata.time(1) <= simtimes(1));
assert(abs(expdata.time(end) - simtimes(end)) < 1e-11);

Evals     = interp1(expdata.time, expdata.U, simtimes, 'linear', 'extrap');
statesExp = cell(numel(output0.states), 1);

for k = 1:numel(output0.states)
    statesExp{k}.time     = simtimes(k);
    statesExp{k}.(ctrl).E = Evals(k);
end

if debug
    % Check that the extracted values are the same as the raw values
    figure; hold on; grid on;
    plot(expdata.time/hour, expdata.U, 'k--');
    plot(getTime(statesExp)/hour, getE(statesExp));
    xlabel('Time / h')
    ylabel('Potential / V')
    title('statesExp')
    drawnow
end

simulatorSetup = SimulationSetup(struct('model'          , output0.model   , ...
                                        'schedule'       , output0.schedule, ...
                                        'initstate'      , output0.initstate, ...
                                        'NonLinearSolver', output0.nls     , ...
                                        'OutputMinisteps', false));

% Setup parameters to be calibrated
HRC = HighRateCalibration(simulatorSetup, 'shortnames', shortnames);
parameters = HRC.getParams();

% Objective function
lsq = @(simsetup, states, varargin) leastSquaresEI(simsetup, states, statesExp, varargin{:});
v = lsq(simulatorSetup, output0.states);
scaling = sum([v{:}]);
objective = @(p, varargin) evalObjectiveBattmo(p, lsq, simulatorSetup, parameters, ...
                                               'objScaling', scaling, varargin{:});

% Compute and classify sensitivities at the initial parameter values.
X0 = getScaledParameterVector(simulatorSetup, parameters);
initialParameterValues = cellfun(@(parameter) parameter.getParameterValue(simulatorSetup), parameters);
sensitivityReport = computeSensitivities(X0, objective, HRC.shortnames, scaling);
senstbl = table(HRC.shortnames(:), abs(sensitivityReport.sensitivities), sensitivityReport.initialGroup(:), ...
                'VariableNames', {'Parameter', 'Sensitivity', 'Initial Group'});
% Sort by sensitivities
[~, sortIdx] = sort(abs(sensitivityReport.sensitivities), 'descend');
senstbl = senstbl(sortIdx, :);
disp(senstbl);
% return

if debug
    % The least squares function evaluated at the experimental values
    % should be zero
    v = lsq(simulatorSetup, statesExp);
    assert(norm([v{:}]) == 0.0);

    % Compare gradients calculated using adjoints and finite
    % difference approximation
    disp('Gradient comparison at initial parameters:');
    compareAdjointAndFiniteDifferenceGradients( ...
        X0, objective, HRC.shortnames, ...
        'PerturbationSize', gradientStepSizes, ...
        'doplot', true);
    %return
end

%% Run optimization

v0 = sensitivityReport.objectiveValue;

callbackfunc = @(history, it) callbackplot(history, it, simulatorSetup, parameters, expdata, ...
                                           'plotEveryIt', 10     , ...
                                           'objScaling' , scaling, ...
                                           'doplot'     , doplot);

gradTol = 1e-3;
objChangeTol = 1e-10;
maxit = 100;
[vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                    'gradTol'         , gradTol     , ...
                                    'objChangeTol'    , objChangeTol, ...
                                    'lineSearchMaxIt' , 10          , ...
                                    'maxInitialUpdate', 0.02        , ...
                                    'maximize'        , false       , ...
                                    'maxit'           , maxit       , ...
                                    'logPlot'         , true        , ...
                                    'callbackfunc'    , callbackfunc, ...
                                    'plotEvolution'   , doplot      , ...
                                    'limitedMemory'   , ~hessian    , ...
                                    'outputHessian'   , hessian);

setupOpt = updateSetupFromScaledParameters(simulatorSetup, parameters, Xopt);

fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));
reasonStr = getReasonStr(history, ...
                         'gradTol'     , gradTol     , ...
                         'objChangeTol', objChangeTol, ...
                         'maxit'       , maxit);
disp(reasonStr);

if debug && numel(history.val) >= 2 && ...
        abs(history.val(end) - history.val(end-1)) < objChangeTol

    % Calculate fd and adjoint gradients at final point
    disp('Gradient comparison at optimized parameters:');
    compareAdjointAndFiniteDifferenceGradients( ...
        Xopt, objective, HRC.shortnames, ...
        'PerturbationSize', gradientStepSizes);
end

% Plot evolution
fig = figure('Position', [100, 100, 560, 560]);
plotParameterEvolution(diaryname, HRC.shortnames(), 'gradTol', gradTol, 'figure', fig);
if dosave
    drawnow
    exportgraphics(fig, '/tmp/parameter-evolution.png', 'resolution', 300)
end

%% Extract parameters

jsonstructHRC = HRC.export(setupOpt);
filename = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
writeStruct(jsonstructHRC, filename);
printer(jsonstructHRC);

%% Run model with calibrated parameters

inputOpt = struct('I'                             , expdata.I                     , ...
                  'totalTime'                     , expdata.time(end)             , ...
                  'numTimesteps'                  , numTimesteps                  , ...
                  'lowRateParams'                 , jsonstructEC                  , ...
                  'highRateParams'                , jsonstructHRC                 , ...
                  'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                  'include_current_collectors'    , true);
outputOpt = runHydra(inputOpt, 'clearSimulation', false);

bgfactorShortnames = {'elyte_bgfactor', 'elyte_bgfactorKappa', 'elyte_bgfactorD'};
bgfactorFields = {'bgfactor', 'bgfactorKappa', 'bgfactorD'};
for k = 1:numel(bgfactorShortnames)
    if any(strcmp(shortnames, bgfactorShortnames{k}))
        fieldname = bgfactorFields{k};
        assert(outputOpt.model.Electrolyte.(fieldname) == ...
               setupOpt.model.Electrolyte.(fieldname));
    end
end

%% Quantify differences

vfinal = lsq(simulatorSetup, outputOpt.states);

expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');
RMSE = l2error(getTime(outputOpt.states), getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);

fprintf('Final least squares values:\n');
fprintf('vopt: %g\n', vopt);
fprintf('Sum of squares: %g\n', sum([vfinal{:}]));
fprintf('RMSE: %g mV\n', RMSE/milli);

if doplot
    % plot differences
    figure; hold on; grid on;
    plot(getTime(outputOpt.states), (getE(outputOpt.states) - expdataUinterp1(getTime(outputOpt.states))).^2, 'displayname', '|E_{sim} - E_{exp}|^2');
    plot(getTime(outputOpt.states), [vfinal{:}], 'displayname', 'vfinal');
end

%% Plot

if doplot
    colors = lines(2);
    fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
    hold on;
    plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2C');
    plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'Initial guess')
    plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'Calibrated');
    xlabel('Time  /  h')
    ylabel('E  /  V')
    legend('location', 'sw')
    axis tight
    ylim([3.45, 4.9])

    if dosave
        exportgraphics(fig, 'high-rate-calibration.png', 'resolution', 300)
    end
end

%% Quantify difference between experiment and calibrated

RMSE = l2error(getTime(outputOpt.states), getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);
fprintf('RMSE after calibration: %g mV\n', RMSE/milli);

%% Check hessian

if hessian

    % history.hess contains the inverse approximate Hessian in scaled coordinates
    invHscaled = history.hess{end};
    issym = @(x) max(max(abs(x-x'))) < 1e-13;
    assert(issym(invHscaled), 'invHscaled is not symmetric');
    symmetrize = @(x) 0.5.*(x + x');
    invHscaled = symmetrize(invHscaled);
    Hscaled = invHscaled \ eye(numel(Xopt));
    assert(issym(Hscaled), 'Hscaled is not symmetric');
    Hscaled = symmetrize(Hscaled);

    % Eigenproblem
    [eigenvecsBFGS, eigenvalsBFGS] = eig(Hscaled, 'vector');
    disp(HRC.shortnames);

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

    % Conditioning and numerical rank
    [U, singularvals, V] = svd(Hscaled, 'econ');
    singularvals = diag(singularvals);

    relativetol = 1e-10;
    ranktol = relativetol * singularvals(1);
    numericalRank = nnz(singularvals > ranktol);
    condno = singularvals(1) / singularvals(end);

    fprintf('Numerical rank: %d/%d\n', numericalRank, numel(Xopt));
    fprintf('SVD condition number: %.3e\n', condno);
    fprintf('Negative eigenvals: %d\n', nnz(eigenvalsBFGS < -ranktol));

    if debug
        % Compare hessian from bfgs with the finite difference of adjoint gradients
        [HfdScaled, HfdReport] = approximateFiniteDifferenceHessian(Xopt, objective, HRC.shortnames, ...
                                                                    'PerturbationSize', hessianStepSizes);

        numberOfStepSizes = numel(hessianStepSizes);
        relativeErrorToBFGS = zeros(numberOfStepSizes, 1);
        maximumAbsoluteError = zeros(numberOfStepSizes, 1);
        fdeigenvals = zeros(numel(Xopt), numberOfStepSizes);

        for stepIndex = 1:numberOfStepSizes
            Hfd = HfdScaled(:, :, stepIndex);
            Hdiff = Hscaled - Hfd;

            relativeErrorToBFGS(stepIndex) = norm(Hdiff, 'fro') ./ max(norm(Hfd, 'fro'), eps);
            maximumAbsoluteError(stepIndex) = max(abs(Hdiff), [], 'all');
            fdeigenvals(:, stepIndex) = sort(eig(Hfd));
        end

        Hcomparison = HfdReport.summary;
        Hcomparison.RelativeErrorToBFGS = relativeErrorToBFGS;
        Hcomparison.MaximumAbsoluteError = maximumAbsoluteError;
        disp('Finite-difference Hessian comparison:');
        disp(Hcomparison);
        fprintf('Gradient evaluations for finite-difference Hessians: %d\n', ...
                HfdReport.numberOfGradientEvaluations);

        stencilNames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(compose('h_%g', hessianStepSizes)));
        stencilTable = cell2table(HfdReport.schemes, ...
                                  'VariableNames', cellstr(stencilNames), ...
                                  'RowNames', HRC.shortnames);
        disp('Finite-difference stencil by parameter:');
        disp(stencilTable);

        eigenvalueComparison = table((1:numel(Xopt))', eigenvalsBFGS, ...
                                     'VariableNames', {'Mode', 'BFGS'});
        for stepIndex = 1:numberOfStepSizes
            varname = matlab.lang.makeValidName(sprintf('FD_h_%g', hessianStepSizes(stepIndex)));
            eigenvalueComparison.(varname) = fdeigenvals(:, stepIndex);
        end
        disp('Hessian eigenvalue comparison:');
        disp(eigenvalueComparison);

        boundtol = 1e-8;
        freeParameters = Xopt > boundtol & Xopt < 1 - boundtol;
        activeParameters = ~freeParameters;

        if any(activeParameters)
            activeParameterTable = table(HRC.shortnames(activeParameters), Xopt(activeParameters), ...
                                         'VariableNames', {'Parameter', 'ScaledValue'});
            disp('Parameters active at a unit-box bound:');
            disp(activeParameterTable);
        end

        fprintf('Free parameters for reduced-Hessian analysis: %d/%d\n', ...
                nnz(freeParameters), numel(Xopt));
        if any(freeParameters)
            for stepIndex = 1:numberOfStepSizes
                Hreduced = HfdScaled(freeParameters, freeParameters, stepIndex);
                reducedEigenvalues = eig(Hreduced);
                fprintf(['FD reduced Hessian h=%g: min eigenvalue=%g, ', ...
                         'max eigenvalue=%g\n'], ...
                        hessianStepSizes(stepIndex), ...
                        min(reducedEigenvalues), max(reducedEigenvalues));
            end
        end
        % keyboard;
    end % end debug

    hessianfdpertsize = 1e-6; % deduced from debug
    [HfdScaled, HfdReport] = approximateFiniteDifferenceHessian(Xopt, objective, HRC.shortnames, ...
                                                                'PerturbationSize', hessianfdpertsize);
    [eigenvecsFD, eigenvalsFD] = eig(HfdScaled, 'vector');

    alleigenvecs = {eigenvecsBFGS, eigenvecsFD};
    alleigenvals = {eigenvalsBFGS, eigenvalsFD};
    casenames = {'BFGS', 'FD'};

    for icase = 1:2

        eigenvecs = alleigenvecs{icase};
        eigenvals = alleigenvals{icase};
        casename = casenames{icase};

        % Plot eigenvectors as a heatmap
        figure;
        imagesc(eigenvecs);
        clim([-1 1]);

        % Diverging colormap: blue -> white -> red
        n = 256;
        c1 = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
        c2 = [ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];

        colormap([c1; c2]);
        colorbar;

        axis equal tight;
        xticks(1:numel(eigenvals));
        modeNames = arrayfun(@(i) sprintf('Mode %d \\lambda_{%d}=%.3g', ...
                                          i, i, eigenvals(i)), ...
                             1:numel(eigenvals), 'UniformOutput', false);
        xticklabels(modeNames);

        yticks(1:numel(eigenvals));
        yticklabels(strrep(HRC.shortnames(), '_', '\_'));

        xlabel(sprintf('%s Hessian eigenmode', casename));
        ylabel('Scaled parameter');
        title(sprintf('Eigenvectors of %s Hessian', casename));
        % set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'tex');

        % Draw grid around the centers
        hold on;
        for k = 0.5:1:(numel(eigenvals)+0.5)
            xline(k, 'k-', 'LineWidth', 0.5);
            yline(k, 'k-', 'LineWidth', 0.5);
        end

        % Add numerical values inside each cell
        for i = 1:numel(eigenvals)
            for j = 1:numel(eigenvals)

                val = eigenvecs(i,j);

                % Choose text color for contrast
                if abs(val) > 0.55
                    txtColor = 'w';
                else
                    txtColor = 'k';
                end

                % Always use two significant digits for clarity
                text(j, i, sprintf('%#.2g', round(val, 2, 'significant')), ...
                     'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', ...
                     'FontWeight', 'bold', ...
                     'FontSize', 13, ...
                     'Color', txtColor);
            end
        end

        if dosave
            drawnow
            exportgraphics(gcf, sprintf('/tmp/hessian-eigenvectors-%s.png', casename), 'resolution', 300)
        end

    end % end icase

end % if hessian

%% Print

disp('Results HRC');
printer(jsonstructHRC);

% Print tortuosities using bgfactor: we have eff cond = bgfactor *
% cond * poro^bman = cond * poro^(bman*lg(bgfactor)) = cond *
% poro^(eff bman). From this we can calculate the tortuosities as tau
% = poro^-bman.

model = outputOpt.model;
if any(strcmp(HRC.shortnames(), 'elyte_bgfactor'))
    rbc = model.(elyte).regionBruggemanCoefficients;
    bgfactor = model.(elyte).bgfactor;
    effbman = @(poro, bman) bman * log(bgfactor) / log(poro);
    effbmen = struct();
    effbmen.(ne) = effbman(1-model.(ne).(co).volumeFraction, rbc.(ne));
    effbmen.(pe) = effbman(1-model.(pe).(co).volumeFraction, rbc.(pe));
    effbmen.(sep) = effbman(model.(sep).porosity, rbc.(sep));
    disp('Pseudo Bruggeman coefficients:');
    printer(effbmen);

    % Convert to tortuosities
    tortuosity = @(vf, bman) vf.^(-bman);
    tau = struct();
    tau.(ne) = tortuosity(1-model.(ne).(co).volumeFraction, effbmen.(ne));
    tau.(pe) = tortuosity(1-model.(pe).(co).volumeFraction, effbmen.(pe));
    tau.(sep) = tortuosity(model.(sep).porosity, effbmen.(sep));
    disp('Derived tortuosities');
    printer(tau);
elseif any(strcmp(HRC.shortnames(), 'elyte_bgfactorKappa')) && any(strcmp(HRC.shortnames(), 'elyte_bgfactorD'))
    rbc = model.(elyte).regionBruggemanCoefficients;
end


% % Postprocess: Report effective electrode conductivities and
% % electrolyte tortuosities
% regionTags = outputOpt.model.(elyte).regionTags;
% vfs = struct();
% vfs.(ne)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 1));
% vfs.(pe)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 2));
% vfs.(sep) = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 3));
% assert(isscalar(vfs.(ne)));
% assert(isscalar(vfs.(pe)));
% assert(isscalar(vfs.(sep)));

% if useRegionBruggemanCoefficients
%     effectiveRegionBruggeman = outputOpt.rbc;
% else
%     bg = outputOpt.model.(elyte).bruggemanCoefficient;
%     assert(isscalar(bg));
%     effectiveRegionBruggeman = struct(ne, bg, pe, bg, sep, bg);
% end

% tortuosityParams = struct();
% tortuosityParams.(elyte).regionBruggemanCoefficients = effectiveRegionBruggeman;
% tau = calculateTortuosityFromBruggeman(vfs, tortuosityParams);
% disp('Effective electrolyte region Bruggeman coefficients:');
% printer(effectiveRegionBruggeman);
% disp('Tortuosities:');
% printer(tau);

effCond = struct(pe, outputOpt.model.(pe).(co).effectiveElectronicConductivity, ...
                 ne, outputOpt.model.(ne).(co).effectiveElectronicConductivity);
disp('Effective electronic conductivities:');
printer(effCond);

fprintf('RMSE %g mV\n', RMSE/milli);
disp(reasonStr)

% Print initial and final vals, plus sensitivities
finalParameterValues = cellfun( ...
    @(parameter) parameter.getParameterValue(setupOpt), parameters);
sensitivitySummary = table( ...
    HRC.shortnames(:), initialParameterValues(:), finalParameterValues(:), ...
    sensitivityReport.absoluteSensitivities(:), sensitivityReport.initialGroup(:), ...
    'VariableNames', ...
    {'Shortname', 'InitialValue', 'FinalValue', 'Sensitivity', 'InitialGroup'});
fprintf('\nInitial sensitivity classification and calibration results:\n');
disp(sensitivitySummary);

diary off;



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
