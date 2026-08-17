function parameter_groups = createClusteredParameterGroups(paramNames, sensitivities, strategy, K)
    % Group parameters using different strategies
    % strategy: 'magnitude', 'physical', 'hybrid_adaptive'

    if nargin < 3
        strategy = 'hybrid_adaptive';
    end
    if nargin < 4
        K = 3; % Default number of groups
    end

    switch strategy
      case 'magnitude'
        parameter_groups = magnitudeBasedGrouping(paramNames, sensitivities);
      case 'physical'
        parameter_groups = physicalDomainGrouping(paramNames, sensitivities);
      case 'hybrid_adaptive'
        parameter_groups = hybridAdaptiveGrouping(paramNames, sensitivities, K);
      otherwise
        parameter_groups = hybridAdaptiveGrouping(paramNames, sensitivities, K);
    end

    fprintf('Grouping strategy: %s\n', strategy);
end

%==========================================================================
% STRATEGY 1: Magnitude-Based Grouping
%==========================================================================
function parameter_groups = magnitudeBasedGrouping(paramNames, sensitivities)
    % Sort by sensitivity values
    [sorted_sens, sort_idx] = sort(abs(sensitivities), 'descend');
    sorted_params = paramNames(sort_idx);

    parameter_groups = {};

    % Simple threshold-based clustering
    if max(sorted_sens) > 0
        % Normalize sensitivities
        norm_sens = sorted_sens / max(sorted_sens);

        % Define closeness thresholds
        high_thresh = 0.2;    % >20% of max sensitivity
        med_thresh = 0.01;    % 1-20% of max sensitivity

        % Group high sensitivity parameters
        high_mask = norm_sens > high_thresh;
        if any(high_mask)
            parameter_groups{end+1} = struct(...
                'name', 'High_Sensitivity', ...
                'parameters', {sorted_params(high_mask)}, ...
                'mean_sensitivity', mean(sorted_sens(high_mask)), ...
                'priority', 1);
        end

        % Group medium sensitivity parameters
        med_mask = norm_sens > med_thresh & norm_sens <= high_thresh;
        if any(med_mask)
            parameter_groups{end+1} = struct(...
                'name', 'Medium_Sensitivity', ...
                'parameters', {sorted_params(med_mask)}, ...
                'mean_sensitivity', mean(sorted_sens(med_mask)), ...
                'priority', 2);
        end

        % Group low sensitivity parameters
        low_mask = norm_sens <= med_thresh;
        if any(low_mask)
            parameter_groups{end+1} = struct(...
                'name', 'Low_Sensitivity', ...
                'parameters', {sorted_params(low_mask)}, ...
                'mean_sensitivity', mean(sorted_sens(low_mask)), ...
                'priority', 3);
        end
    else
        % Fallback: equal groups
        n_params = length(sorted_params);
        group_size = ceil(n_params / 3);
        for i = 1:3
            start_idx = (i-1) * group_size + 1;
            end_idx = min(i * group_size, n_params);
            if start_idx <= n_params
                parameter_groups{end+1} = struct(...
                    'name', sprintf('Group_%d', i), ...
                    'parameters', {sorted_params(start_idx:end_idx)}, ...
                    'mean_sensitivity', 0, ...
                    'priority', i);
            end
        end
    end

    displayResults(parameter_groups, 'Magnitude-Based');
end

%==========================================================================
% STRATEGY 2: Physical Domain Grouping
%==========================================================================
function parameter_groups = physicalDomainGrouping(paramNames, sensitivities)
    parameter_groups = {};

    % Define physical domains based on battery model physics
    domains = {
        {'Electrode_Kinetics', {'ne_vsa', 'pe_vsa'}};
        {'Electrode_Transport', {'ne_bg', 'pe_bg'}};
        {'Electrolyte_Transport', {'elyte_bg'}};
        {'Solid_Diffusion', {'ne_D', 'pe_D'}};
        {'Conductivity', {'ne_kappa', 'pe_kappa'}}
              };

    % We should not calibrate conductivity
    assert(~ismember('kappa', paramNames), 'Do not calibrate conductivity parameters');

    abs_sensitivities = abs(sensitivities);

    % Create groups for each domain with available parameters
    priority = 1;
    for i = 1:length(domains)
        domain_name = domains{i}{1};
        domain_params = domains{i}{2};
        available_params = intersect(paramNames, domain_params);

        if ~isempty(available_params)
            domain_sens = abs_sensitivities(ismember(paramNames, available_params));

            parameter_groups{end+1} = struct(...
                'name', domain_name, ...
                'parameters', {available_params}, ...
                'mean_sensitivity', mean(domain_sens), ...
                'priority', priority);
            priority = priority + 1;
        end
    end

    % Handle any ungrouped parameters
    all_grouped = {};
    for i = 1:length(parameter_groups)
        all_grouped = [all_grouped, parameter_groups{i}.parameters];
    end
    ungrouped = setdiff(paramNames, all_grouped);

    if ~isempty(ungrouped)
        ungrouped_sens = abs_sensitivities(ismember(paramNames, ungrouped));
        parameter_groups{end+1} = struct(...
            'name', 'Ungrouped', ...
            'parameters', {ungrouped}, ...
            'mean_sensitivity', mean(ungrouped_sens), ...
            'priority', priority);
    end

    % Sort by mean sensitivity (highest first)
    if ~isempty(parameter_groups)
        mean_sens = cellfun(@(group) group.mean_sensitivity, parameter_groups);
        [~, sort_idx] = sort(mean_sens, 'descend');
        parameter_groups = parameter_groups(sort_idx);

        % Update priorities
        for i = 1:length(parameter_groups)
            parameter_groups{i}.priority = i;
        end
    end

    displayResults(parameter_groups, 'Physical Domain');
end

%==========================================================================
% STRATEGY 3: Hybrid Adaptive Grouping (RECOMMENDED)
%==========================================================================
function parameter_groups = hybridAdaptiveGrouping(paramNames, sensitivities, K)
    abs_sensitivities = abs(sensitivities);

    % Adaptive threshold selection
    max_sens = max(abs_sensitivities);
    min_sens = min(abs_sensitivities);
    sensitivity_ratio = max_sens / (min_sens + eps);

    if sensitivity_ratio > 1000
        % Large dynamic range - use logarithmic percentiles
        log_sens = log10(abs_sensitivities + eps);
        high_threshold = 10^prctile(log_sens, 70);
        medium_threshold = 10^prctile(log_sens, 30);
    else
        % Moderate range - use linear percentiles
        high_threshold = prctile(abs_sensitivities, 70);
        medium_threshold = prctile(abs_sensitivities, 30);
    end

    % Create groups
    high_mask = abs_sensitivities >= high_threshold;
    medium_mask = abs_sensitivities >= medium_threshold & abs_sensitivities < high_threshold;
    low_mask = abs_sensitivities < medium_threshold;

    parameter_groups = {};

    % High sensitivity group
    high_params = paramNames(high_mask);
    if ~isempty(high_params)
        high_sens = sensitivities(high_mask);
        [~, high_sort] = sort(abs(high_sens), 'descend');

        parameter_groups{end+1} = struct(...
            'name', 'High_Sensitivity', ...
            'parameters', {high_params(high_sort)}, ...
            'mean_sensitivity', mean(abs(high_sens)), ...
            'priority', 1);
    end

    % Medium sensitivity group
    medium_params = paramNames(medium_mask);
    if ~isempty(medium_params)
        medium_sens = sensitivities(medium_mask);
        [~, med_sort] = sort(abs(medium_sens), 'descend');

        parameter_groups{end+1} = struct(...
            'name', 'Medium_Sensitivity', ...
            'parameters', {medium_params(med_sort)}, ...
            'mean_sensitivity', mean(abs(medium_sens)), ...
            'priority', 2);
    end

    % Low sensitivity group
    low_params = paramNames(low_mask);
    if ~isempty(low_params)
        low_sens = sensitivities(low_mask);

        parameter_groups{end+1} = struct(...
            'name', 'Low_Sensitivity', ...
            'parameters', {low_params}, ...
            'mean_sensitivity', mean(abs(low_sens)), ...
            'priority', 3);
    end

    displayResults(parameter_groups, 'Hybrid Adaptive');
end

%==========================================================================
% Common display function
%==========================================================================
function displayResults(parameter_groups, strategy_name)
    fprintf('\n%s Grouping Results:\n', strategy_name);

    for i = 1:length(parameter_groups)
        group = parameter_groups{i};
        fprintf('\n%s (Priority %d):\n', group.name, group.priority);
        fprintf('  Mean sensitivity: %.2e\n', group.mean_sensitivity);

        % Safe parameter display for any data type
        if isfield(group, 'parameters')
            if iscell(group.parameters)
                fprintf('  Parameters (%d): ', numel(group.parameters));
                for j = 1:numel(group.parameters)
                    if j > 1, fprintf(', '); end
                    fprintf('%s', group.parameters{j});
                end
            elseif ischar(group.parameters)
                fprintf('  Parameter: %s', group.parameters);
            else
                fprintf('  Parameters: [array of %d elements]', numel(group.parameters));
            end
        else
            fprintf('  No parameters field');
        end
        fprintf('\n');
    end

    fprintf('\nRecommended Estimation Sequence:\n');
    for i = 1:length(parameter_groups)
        fprintf('  Step %d: %s\n', parameter_groups{i}.priority, parameter_groups{i}.name);
    end
end
