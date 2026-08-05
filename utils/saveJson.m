function saveJson(filename, data)

    jsonstr = jsonencode(data, 'PrettyPrint', true);
    fid = fopen(filename, 'w');
    assert(fid > 0, 'Could not open %s for writing', filename);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s', jsonstr);

end
