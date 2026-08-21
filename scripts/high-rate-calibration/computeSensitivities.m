function report = computeSensitivities(scaledParameters, objective, ...
                                       shortnames, objectiveScaling)
% Compute and classify the initial (high-rate) calibration sensitivities.

    assert(isa(objective, 'function_handle'), ...
           'objective must be a function handle.');
    assert(isscalar(objectiveScaling) && isfinite(objectiveScaling), ...
           'objectiveScaling must be a finite scalar.');

    shortnames = reshape(shortnames, [], 1);
    [objectiveValue, scaledGradient] = objective(scaledParameters);
    sensitivities = reshape(scaledGradient .* objectiveScaling, [], 1);

    assert(numel(shortnames) == numel(sensitivities), ...
           'Expected one sensitivity for each parameter shortname.');
    assert(all(isfinite(sensitivities)), ...
           'All parameter sensitivities must be finite.');

    initialGroup = classifySensitivities(sensitivities);

    fprintf('\n=== Initial parameter sensitivities ===\n');
    for index = 1:numel(shortnames)
        level = erase(initialGroup{index}, '_Sensitivity');
        fprintf('  %-20s % .3e  %s sensitivity\n', ...
                shortnames{index}, abs(sensitivities(index)), level);
    end

    report = struct( ...
        'objectiveValue', objectiveValue, ...
        'scaledGradient', scaledGradient, ...
        'sensitivities', sensitivities, ...
        'initialGroup', {initialGroup});

end


function initialGroup = classifySensitivities(sensitivities)
% Hybrid-adaptive grouping copied from the sequential calibration workflow.

    absoluteSensitivities = abs(sensitivities);
    initialGroup = repmat({''}, numel(absoluteSensitivities), 1);

    if all(absoluteSensitivities == 0)
        initialGroup(:) = {'Low_Sensitivity'};
        return
    end

    maxSensitivity = max(absoluteSensitivities);
    minSensitivity = min(absoluteSensitivities);
    sensitivityRatio = maxSensitivity / (minSensitivity + eps);

    if sensitivityRatio > 1000
        logSensitivities = log10(absoluteSensitivities + eps);
        highThreshold = 10^prctile(logSensitivities, 70);
        mediumThreshold = 10^prctile(logSensitivities, 30);
    else
        highThreshold = prctile(absoluteSensitivities, 70);
        mediumThreshold = prctile(absoluteSensitivities, 30);
    end

    high = absoluteSensitivities >= highThreshold;
    medium = absoluteSensitivities >= mediumThreshold & ...
             absoluteSensitivities < highThreshold;
    low = absoluteSensitivities < mediumThreshold;

    initialGroup(high) = {'High_Sensitivity'};
    initialGroup(medium) = {'Medium_Sensitivity'};
    initialGroup(low) = {'Low_Sensitivity'};

    assert(all(~cellfun(@isempty, initialGroup)), ...
           'Every parameter must receive a sensitivity classification.');

end
