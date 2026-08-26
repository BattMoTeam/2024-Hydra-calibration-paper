%% Plot optimizer convergence and parameter evolution from the calibration log

clearvars
close all

scriptDirectory = fileparts(mfilename('fullpath'));
logFilename = fullfile(scriptDirectory, ...
    '_diary-runHighRateCalibration-20260821-103811.txt');

assert(isfile(logFilename), 'Calibration log not found: %s', logFilename);
logText = fileread(logFilename);

numberPattern = '[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?';

% pgrad is the infinity norm of the projected objective gradient. It is a
% scalar convergence diagnostic, not an individual parameter sensitivity.
pgradPattern = ['It:\s*(\d+)\s*\|[^\r\n]*?pgrad:\s*(' ...
    numberPattern ')'];
pgradTokens = regexp(logText, pgradPattern, 'tokens');
assert(~isempty(pgradTokens), 'No pgrad records found in %s.', logFilename);

pgradIterations = cellfun(@(entry) str2double(entry{1}), pgradTokens).';
pgradValues = cellfun(@(entry) str2double(entry{2}), pgradTokens).';
assert(all(isfinite(pgradIterations)) && all(isfinite(pgradValues)), ...
    'The log contains invalid iteration or pgrad values.');
assert(all(diff(pgradIterations) > 0), ...
    'Expected pgrad iterations to be strictly increasing.');

% The scaled objective functional is logged at every optimizer iteration.
objectivePattern = ['It:\s*(\d+)\s*\|\s*val:\s*(' ...
    numberPattern ')'];
objectiveTokens = regexp(logText, objectivePattern, 'tokens');
assert(~isempty(objectiveTokens), ...
    'No objective-functional records found in %s.', logFilename);

objectiveIterations = cellfun( ...
    @(entry) str2double(entry{1}), objectiveTokens).';
objectiveValues = cellfun( ...
    @(entry) str2double(entry{2}), objectiveTokens).';
assert(all(isfinite(objectiveIterations)) && ...
    all(isfinite(objectiveValues)), ...
    'The log contains invalid iteration or objective values.');
assert(all(objectiveValues > 0) && ...
    isequal(objectiveIterations, pgradIterations), ...
    'Expected a positive objective value for every pgrad iteration.');

% callbackplot logs the physical parameter values after every completed
% optimizer iteration. Iteration zero is reconstructed from its repeated
% "initial values" record.
parameterPattern = [ ...
    'callbackplot it=(\d+)\r?\n' ...
    'vad[^\r\n]*\r?\n' ...
    'u[^\r\n]*\r?\n' ...
    'initial values\s+([^\r\n]+)\r?\n' ...
    'vals\s+([^\r\n]+)'];
parameterTokens = regexp(logText, parameterPattern, 'tokens');
assert(~isempty(parameterTokens), ...
    'No callback parameter records found in %s.', logFilename);

initialValues = sscanf(parameterTokens{1}{2}, '%f').';
numberOfParameters = numel(initialValues);
callbackIterations = zeros(numel(parameterTokens), 1);
callbackValues = zeros(numel(parameterTokens), numberOfParameters);

for recordIndex = 1:numel(parameterTokens)
    callbackIterations(recordIndex) = str2double(parameterTokens{recordIndex}{1});
    values = sscanf(parameterTokens{recordIndex}{3}, '%f').';
    assert(numel(values) == numberOfParameters, ...
        'Iteration %d has %d parameter values; expected %d.', ...
        callbackIterations(recordIndex), numel(values), numberOfParameters);
    callbackValues(recordIndex, :) = values;
end

assert(all(isfinite(initialValues)) && all(isfinite(callbackValues), 'all'), ...
    'The log contains invalid parameter values.');
assert(all(initialValues ~= 0), ...
    'Cannot normalize parameters whose initial value is zero.');
assert(all(diff(callbackIterations) > 0), ...
    'Expected callback iterations to be strictly increasing.');

parameterIterations = [0; callbackIterations];
parameterValues = [initialValues; callbackValues];
relativeParameterValues = parameterValues ./ initialValues;

parameterNames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', ...
    'elyte_bgfactor'};
assert(numberOfParameters == numel(parameterNames), ...
    'Found %d parameters, but %d names are configured.', ...
    numberOfParameters, numel(parameterNames));

colors = lines(numberOfParameters);
figureHandle = figure('Color', 'w', ...
    'Name', 'High-rate calibration parameter evolution');
layout = tiledlayout(figureHandle, 3, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

gradientAxes = nexttile(layout);
semilogy(gradientAxes, pgradIterations, pgradValues, 'o-', ...
    'Color', [0.15, 0.15, 0.15], 'MarkerFaceColor', [0.15, 0.15, 0.15], ...
    'MarkerSize', 4, 'DisplayName', 'Projected gradient norm');
grid(gradientAxes, 'on');
ylabel(gradientAxes, 'Projected gradient norm');
legend(gradientAxes, 'Location', 'northeast');

objectiveAxes = nexttile(layout);
semilogy(objectiveAxes, objectiveIterations, objectiveValues, 's-', ...
    'Color', [0.00, 0.45, 0.74], ...
    'MarkerFaceColor', [0.00, 0.45, 0.74], 'MarkerSize', 5, ...
    'DisplayName', 'Objective functional');
grid(objectiveAxes, 'on');
ylabel(objectiveAxes, 'Objective functional');
legend(objectiveAxes, 'Location', 'northeast');

parameterAxes = nexttile(layout);
hold(parameterAxes, 'on');
for parameterIndex = 1:numberOfParameters
    plot(parameterAxes, parameterIterations, ...
        relativeParameterValues(:, parameterIndex), 'o-', ...
        'Color', colors(parameterIndex, :), 'MarkerSize', 4, ...
        'DisplayName', parameterNames{parameterIndex});
end
grid(parameterAxes, 'on');
xlabel(parameterAxes, 'Optimization iteration');
ylabel(parameterAxes, 'Parameter value / initial value');
legend(parameterAxes, 'Location', 'eastoutside', 'Interpreter', 'none');

title(layout, 'High-rate calibration evolution');
linkaxes([gradientAxes, objectiveAxes, parameterAxes], 'x');
xlim(gradientAxes, [min(pgradIterations), max(pgradIterations)]);
