function saveFigureSet(fig, outbase)

    [folder, ~, ~] = fileparts(outbase);
    ensureFolder(folder);
    savefig(fig, [outbase, '.fig']);
    exportgraphics(fig, [outbase, '.png'], 'Resolution', 300);

end
