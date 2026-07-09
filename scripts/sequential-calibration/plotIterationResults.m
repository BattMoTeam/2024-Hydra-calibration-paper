function plotIterationResults(iteration, output_ref, current_output, expdata, gradient_overall_actual, paramNames, optimization_history, varargin)

    opt = struct('tag', '');
    opt = merge_options(opt, varargin{:});

    fig = figure(100 + iteration);
    clf(fig);
    set(fig, 'Position', [100, 100, 1400, 1000]);

    ctrl = 'Control';
    getTime = @(states) cellfun(@(s) s.time, states);
    getE = @(states) cellfun(@(s) s.(ctrl).E, states);
    lineWidth = 2;
    plotColors = lines(2);
    calibratedColor = plotColors(1, :);
    initialColor = plotColors(2, :);
    experimentalColor = 'k';

    % 1. Voltage comparison
    subplot(2, 3, 1);
    plot(expdata.time/hour, expdata.E, '--', 'Color', experimentalColor, ...
         'LineWidth', lineWidth, 'DisplayName', 'Experimental');
    hold on;
    plot(getTime(output_ref.states)/hour, getE(output_ref.states), '-', ...
         'Color', initialColor, 'LineWidth', lineWidth, 'DisplayName', 'Initial');
    plot(getTime(current_output.states)/hour, getE(current_output.states), '-', ...
         'Color', calibratedColor, 'LineWidth', lineWidth, 'DisplayName', 'Calibrated');
    xlabel('Time  /  h');
    ylabel('Voltage  /  V');
    title('Voltage Comparison', 'Interpreter', 'none');
    legend('Location', 'best');
    grid on;

    % 2. Voltage error
    subplot(2, 3, 2);
    E_sim_initial = interp1(getTime(output_ref.states), getE(output_ref.states), expdata.time, 'linear', 'extrap');
    E_sim_optimized = interp1(getTime(current_output.states), getE(current_output.states), expdata.time, 'linear', 'extrap');

    error_initial = (E_sim_initial - expdata.E) * 1000;
    error_optimized = (E_sim_optimized - expdata.E) * 1000;

    plot(expdata.time/hour, error_initial, '-', ...
         'Color', initialColor, 'LineWidth', lineWidth, 'DisplayName', 'Initial');
    hold on;
    plot(expdata.time/hour, error_optimized, '-', ...
         'Color', calibratedColor, 'LineWidth', lineWidth, 'DisplayName', 'Calibrated');
    xlabel('Time  /  h');
    ylabel('Error  /  mV');
    title('Voltage Error', 'Interpreter', 'none');
    legend('Location', 'best');
    grid on;

    % 3. Parameter sensitivities
    subplot(2, 3, 3);
    sensitivities = abs(gradient_overall_actual(:));
    [sensitivities_sorted, sort_idx] = sort(sensitivities, 'descend');
    paramNames_sorted = paramNames(sort_idx);
    bar(sensitivities_sorted, 'FaceColor', calibratedColor);
    set(gca, 'XTick', 1:length(paramNames_sorted), ...
             'XTickLabel', strrep(paramNames_sorted, '_', '\_'), ...
             'XTickLabelRotation', 45);
    ylabel('Sensitivity');
    title('Parameter Sensitivities', 'Interpreter', 'none');
    grid on;
    hold on;
    group_ids = getSortedGroupIds(paramNames_sorted, optimization_history);
    addGroupSeparators(gca, group_ids);

    % 4. Objective function history
    subplot(2, 3, 4);
    if ~isempty(optimization_history)
        group_start = 1;
        for i = 1:length(optimization_history)
            if isfield(optimization_history(i), 'iteration_history')
                objective = optimization_history(i).iteration_history;
                if ~isempty(objective)
                    objective = objective ./ objective(1);
                    x = group_start:(group_start + numel(objective) - 1);
                    plot(x, objective, '-o', 'Color', calibratedColor, ...
                         'LineWidth', lineWidth, 'MarkerFaceColor', calibratedColor);
                    hold on;
                    if i > 1
                        xline(group_start - 1, 'k--', 'LineWidth', 1);
                    end
                    yl = ylim(gca);
                    text(mean(x), labelYPosition(gca, 1.5), sprintf('G%d', i), ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'top', ...
                         'FontWeight', 'bold');
                    ylim(yl);
                    group_start = group_start + numel(objective) + 1;
                end
            end
        end
        if group_start > 1
            xlabel('Function Evaluations');
            ylabel('Normalized Objective Value');
            title('Optimization Progress', 'Interpreter', 'none');
            grid on;
        end
    end

    % 5. Group improvements
    subplot(2, 3, 5);
    if ~isempty(optimization_history)
        improvements = [];
        group_labels = {};
        for i = 1:length(optimization_history)
            if isfield(optimization_history(i), 'improvement')
                improvements(end+1) = optimization_history(i).improvement;
                group_labels{end+1} = sprintf('G%d', i);
            end
        end
        if ~isempty(improvements)
            bar(improvements, 'FaceColor', calibratedColor);
            set(gca, 'XTickLabel', group_labels);
            ylabel('Improvement  /  %');
            title('Group Improvements', 'Interpreter', 'none');
            grid on;
        end
    end

    % 6. L2 error comparison
    subplot(2, 3, 6);
    initial_wL2 = l2error(getTime(output_ref.states), getE(output_ref.states), ...
                          expdata.time, expdata.E, 'extrap', true);
    final_wL2 = l2error(getTime(current_output.states), getE(current_output.states), ...
                        expdata.time, expdata.E, 'extrap', true);

    b = bar([initial_wL2/milli, final_wL2/milli]);
    b.FaceColor = 'flat';
    b.CData = [initialColor; calibratedColor];
    set(gca, 'XTickLabel', {'Initial', 'Final'});
    ylabel('RMSE  /  mV');
    title('RMSE Comparison', 'Interpreter', 'none');
    grid on;

    % Add L2 values as text
    text(1, initial_wL2/milli, sprintf('%2.2f', initial_wL2/milli), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    text(2, final_wL2/milli, sprintf('%2.2f', final_wL2/milli), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    sgtitle(sprintf('Iteration %d Results', iteration), 'FontSize', 16, 'FontWeight', 'bold');

    saveas(fig, sprintf('%siteration_%d_results.png', opt.tag, iteration));
end

function group_ids = getSortedGroupIds(paramNames_sorted, optimization_history)
    group_ids = nan(size(paramNames_sorted));
    if isempty(optimization_history)
        return
    end
    for igroup = 1:numel(optimization_history)
        if ~isfield(optimization_history(igroup), 'parameters')
            continue
        end
        group_params = optimization_history(igroup).parameters;
        for iparam = 1:numel(paramNames_sorted)
            if any(strcmp(paramNames_sorted{iparam}, group_params))
                group_ids(iparam) = igroup;
            end
        end
    end
end

function addGroupSeparators(ax, group_ids)
    if isempty(group_ids) || all(isnan(group_ids))
        return
    end

    yl = ylim(ax);
    run_start = 1;
    for i = 2:(numel(group_ids) + 1)
        if i <= numel(group_ids)
            same_group = isequaln(group_ids(i), group_ids(i - 1));
        else
            same_group = false;
        end

        if same_group
            continue
        end

        run_end = i - 1;
        group_id = group_ids(run_start);
        if ~isnan(group_id)
            xline(ax, run_start - 0.5, 'k--', 'LineWidth', 1);
            xline(ax, run_end + 0.5, 'k--', 'LineWidth', 1);
            text(ax, mean([run_start, run_end]), labelYPosition(ax, 1.5), sprintf('G%d', group_id), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'top', ...
                 'FontWeight', 'bold');
        end
        run_start = i;
    end
end

function y = labelYPosition(ax, characterHeights)
    yl = ylim(ax);
    fontSizePoints = ax.FontSize;
    axesHeightPoints = ax.Position(4) * ax.Parent.Position(4);
    if axesHeightPoints <= 0
        y = yl(2);
        return
    end
    labelOffset = characterHeights * fontSizePoints / axesHeightPoints * diff(yl);
    y = yl(2) - labelOffset;
end
