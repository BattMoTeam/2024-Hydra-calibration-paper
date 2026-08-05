function values = stackStateVectors(states, getter)

    vectors = cellfun(@(s) getter(s), states, 'UniformOutput', false);
    values = horzcat(vectors{:});

end
