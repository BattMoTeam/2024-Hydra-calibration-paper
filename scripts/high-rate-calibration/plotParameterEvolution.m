function plotParameterEvolution(logFilename, shortnames, varargin)
% Plot convergence and parameter values recorded by callbackplot.

    opt = struct('gradTol', [], ...
                 'figure', []);
    opt = merge_options(opt, varargin{:});

    assert(isfile(logFilename), 'Calibration log not found: %s', logFilename);
    txt = fileread(logFilename);

    number = '[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?';
    getcol = @(tokens, col) cellfun(@(x) str2double(x{col}), tokens).';

    tokens = regexp(txt, ['It:\s*(\d+)\s*\|[^\r\n]*?pgrad:\s*(' number ')'], 'tokens');
    assert(~isempty(tokens), 'No pgrad records found in %s.', logFilename);
    iterations = getcol(tokens, 1);
    pgrad = getcol(tokens, 2);

    tokens = regexp(txt, ['It:\s*(\d+)\s*\|\s*val:\s*(' number ')'], 'tokens');
    assert(~isempty(tokens), 'No objective records found in %s.', logFilename);
    objectiveIterations = getcol(tokens, 1);
    objective = getcol(tokens, 2);

    assert(isequal(objectiveIterations, iterations), ...
        'Objective and gradient records do not use the same iterations.');
    assert(all(objective > 0), 'Objective values must be positive for semilogy.');

    pattern = ['callbackplot it=(\d+)\r?\n' ...
               'vad[^\r\n]*\r?\n' ...
               'u[^\r\n]*\r?\n' ...
               'initial values\s+([^\r\n]+)\r?\n' ...
               'vals\s+([^\r\n]+)'];
    tokens = regexp(txt, pattern, 'tokens');
    assert(~isempty(tokens), 'No callback records found in %s.', logFilename);

    initial = sscanf(tokens{1}{2}, '%f').';
    nparam = numel(initial);
    callbackIterations = zeros(numel(tokens), 1);
    values = zeros(numel(tokens), nparam);

    for i = 1:numel(tokens)
        callbackIterations(i) = str2double(tokens{i}{1});
        row = sscanf(tokens{i}{3}, '%f').';
        assert(numel(row) == nparam, ...
            'Iteration %d has the wrong number of parameters.', callbackIterations(i));
        values(i, :) = row;
    end

    assert(all(initial ~= 0), 'Cannot normalize a parameter with initial value zero.');
    assert(nparam == numel(shortnames), ...
        'Found %d parameters but received %d names.', nparam, numel(shortnames));

    parameterIterations = [0; callbackIterations];
    relativeValues = [initial; values] ./ initial;

    if isempty(opt.figure)
        fig = figure;
    else
        fig = opt.figure;
    end

    layout = tiledlayout(fig, 3, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    markerSize = 18;

    ax(1) = nexttile(layout);
    semilogy(ax(1), iterations, pgrad, 'k.-', 'MarkerSize', markerSize);
    grid(ax(1), 'on');
    ylabel(ax(1), 'Projected gradient norm');
    if ~isempty(opt.gradTol)
        yline(ax(1), opt.gradTol, 'k--');
    end

    ax(2) = nexttile(layout);
    semilogy(ax(2), objectiveIterations, objective, 'k.-', ...
        'MarkerSize', markerSize);
    grid(ax(2), 'on');
    ylabel(ax(2), 'Objective functional');

    ax(3) = nexttile(layout);
    hold(ax(3), 'on');
    colors = lines(nparam);
    for i = 1:nparam
        plot(ax(3), parameterIterations, relativeValues(:, i), '.-', ...
            'Color', colors(i, :), 'MarkerSize', markerSize, ...
            'DisplayName', shortnames{i});
    end
    grid(ax(3), 'on');
    xlabel(ax(3), 'Optimization iteration');
    ylabel(ax(3), 'Normalized param. value');
    legend(ax(3), 'Location', 'northwest', 'Interpreter', 'none');

    title(layout, 'High-rate calibration evolution');
    linkaxes(ax, 'x');
    xlim(ax(1), [min(iterations), max(iterations)]);

end
