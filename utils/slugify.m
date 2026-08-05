function slug = slugify(str)

    slug = lower(strrep(str, ' ', '-'));
    slug = regexprep(slug, '[^a-z0-9\-]', '');

end
