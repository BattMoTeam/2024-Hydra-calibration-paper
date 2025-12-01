function fjv = jsonDiff(json1input, json2input, name1, name2)

    if ~isstr(name1)
        name1 = sprintf('%s', name1);
    end

    if ~isstr(name2)
        name2 = sprintf('%s', name2);
    end

    if ischar(json1input)
        json1 = jsondecode(fileread(json1input));
    else
        json1 = json1input;
    end

    if ischar(json2input)
        json2 = jsondecode(fileread(json2input));
    else
        json2 = json2input;
    end

    fjv = compareJson(json1, json2, name1, name2);
    fjv.print('filter', {'comparison', 'different'});

end
