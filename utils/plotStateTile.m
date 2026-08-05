function plotStateTile(ax, time_h, x_um, values, plotTitle)

    axes(ax);
    imagesc(time_h, x_um, values);
    axis xy
    grid off
    xlabel('Time / h')
    ylabel('x / \mum')
    title(plotTitle)
    colorbar

end
