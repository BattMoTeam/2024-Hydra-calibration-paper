function ensureFolder(folder)

    if ~isfolder(folder)
        mkdir(folder);
    end

end
