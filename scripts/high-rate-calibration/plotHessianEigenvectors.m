function fig = plotHessianEigenvectors(H, shortnames, casename, varargin)
% Plot the eigenvectors of a Hessian as a heatmap.

    opt = struct('dosave', false);
    opt = merge_options(opt, varargin{:});

    H = 0.5.*(H + H');
    [eigenvecs, eigenvals] = eig(H, 'vector');

    for k = 1:numel(eigenvals)
        [~, maxIdx] = max(abs(eigenvecs(:, k)));
        if eigenvecs(maxIdx, k) < 0
            eigenvecs(:, k) = -eigenvecs(:, k);
        end
    end

    fig = figure;
    imagesc(eigenvecs);
    clim([-1 1]);

    n = 256;
    c1 = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
    c2 = [ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];
    colormap([c1; c2]);
    colorbar;

    axis equal tight;
    xticks(1:numel(eigenvals));
    modeNames = arrayfun(@(idx) sprintf('Mode %d \\lambda_{%d}=%.3g', ...
                                      idx, idx, eigenvals(idx)), ...
                         1:numel(eigenvals), 'UniformOutput', false);
    xticklabels(modeNames);
    yticks(1:numel(eigenvals));
    yticklabels(strrep(shortnames, '_', '\_'));

    xlabel(sprintf('%s Hessian eigenmode', casename));
    ylabel('Scaled parameter');
    title(sprintf('Eigenvectors of %s Hessian', casename));

    hold on;
    for k = 0.5:1:(numel(eigenvals)+0.5)
        xline(k, 'k-', 'LineWidth', 0.5);
        yline(k, 'k-', 'LineWidth', 0.5);
    end

    for i = 1:numel(eigenvals)
        for j = 1:numel(eigenvals)
            val = eigenvecs(i,j);
            if abs(val) > 0.55
                txtColor = 'w';
            else
                txtColor = 'k';
            end

            text(j, i, sprintf('%#.2g', round(val, 2, 'significant')), ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'FontWeight', 'bold', ...
                 'FontSize', 13, ...
                 'Color', txtColor);
        end
    end

    if opt.dosave
        drawnow
        exportgraphics(fig, sprintf('/tmp/hessian-eigenvectors-%s.png', casename), ...
                       'resolution', 300)
    end

end
