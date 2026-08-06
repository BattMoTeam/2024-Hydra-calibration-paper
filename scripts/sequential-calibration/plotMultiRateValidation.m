function plotMultiRateValidation(expdata_all, lowRateParams, highRateParams, iteration, varargin)

    opt = struct('tag'                           , '', ...
                 'geometry'                      , '1d', ...
                 'include_current_collectors'    , false, ...
                 'useRegionBruggemanCoefficients', false);
    opt = merge_options(opt, varargin{:});

% Plot validation across multiple discharge rates
    fig = figure(200 + iteration);
    clf(fig);
    set(fig, 'Position', [100, 100, 1200, 800]);
    ax = axes(fig);

    fig2 = figure(300 + iteration);
    clf(fig2);
    set(fig2, 'Position', [150, 150, 1200, 800]);
    ax2 = axes(fig2);
    hold(ax2, 'on');
    grid(ax2, 'on');
    legend(ax2, 'Location', 'sw');
    xlabel(ax2, 'Capacity  /  Ah');
    ylabel(ax2, 'Voltage  /  V');

    colors = lines(numel(expdata_all));

    ncols = ceil(sqrt(numel(expdata_all)));
    nrows = ceil(numel(expdata_all) / ncols);

    % Select 4 different rates for validation
    %validation_rates = [1, 2, 3, 4]; % Adjust indices based on your expdata_all
    validation_rates = 1:numel(expdata_all);
    % if length(validation_rates) > length(expdata_all)
    %     validation_rates = 1:min(4, length(expdata_all));
    % end

    L2_errors = zeros(1, length(validation_rates));
    wL2_errors = zeros(1, length(validation_rates));

    ctrl = 'Control';
    getTime = @(states) cellfun(@(s) s.time, states);
    getE = @(states) cellfun(@(s) s.(ctrl).E, states);
    getI = @(states) cellfun(@(s) s.(ctrl).I, states);

    for i = 1:length(validation_rates)
        rate_idx = validation_rates(i);
        current_expdata = expdata_all{rate_idx};

        % Create input for this rate
        input_val = struct('I', current_expdata.I, ...
                           'totalTime', current_expdata.time(end), ...
                           'geometry', opt.geometry, ...
                           'include_current_collectors', opt.include_current_collectors, ...
                           'lowRateParams', lowRateParams, ...
                           'highRateParams', highRateParams, ...
                           'useRegionBruggemanCoefficients', opt.useRegionBruggemanCoefficients);

        % Run simulation with optimized parameters
        output_val = runHydra(input_val, 'clearSimulation', false);

        % Calculate L2 error
        wL2_errors(i) = l2error(getTime(output_val.states), getE(output_val.states), ...
                                current_expdata.time, current_expdata.E, ...
                                'extrap', true);

        % Plot
        figure(fig);
        subplot(nrows, ncols, i);
        plot(current_expdata.time/hour, current_expdata.E, 'k-', 'LineWidth', 2, 'DisplayName', 'Experimental');
        hold on;
        plot(getTime(output_val.states)/hour, getE(output_val.states), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Simulated');

        xlabel('Time  /  h');
        ylabel('Voltage  /  V');
        str = sprintf('Rate %.3fC (RMSE = %2.2f mV)', current_expdata.rate, wL2_errors(i)/milli);
        title(str, 'FontSize', 12);
        legend('Location', 'best');
        grid on;

        plot(ax2, current_expdata.Q/hour, current_expdata.E, '--', 'LineWidth', 2, 'DisplayName', sprintf('Exp %.3fC', current_expdata.rate), 'Color', colors(i, :));
        Qsim = cumtrapz(getTime(output_val.states), getI(output_val.states));
        plot(ax2, Qsim/hour, getE(output_val.states), '-', 'LineWidth', 1.5, 'DisplayName', strrep(str, 'Rate', 'Sim'), 'Color', colors(i, :));

    end

    sgtitle(sprintf('Multi-Rate Validation - Iteration %d', iteration), 'FontSize', 14, 'FontWeight', 'bold');

    % Save figure
    saveas(fig, sprintf('%smultirate_validation_iter_%d.png', opt.tag, iteration));

    % Display L2 errors
    fprintf('\nMulti-Rate Validation RMSE:\n');
    for i = 1:length(validation_rates)
        rate_idx = validation_rates(i);
        fprintf('  Rate %.3fC: RMSE = %2.2f mV\n', expdata_all{rate_idx}.rate, wL2_errors(i)/milli);
    end
    fprintf('  Average RMSE: %2.2f mV\n', mean(wL2_errors)/milli);
end
