function volumeFractions(paramobj)
% Set up volume fractions from mass fractions

    ne = 'NegativeElectrode';
    pe = 'PositiveElectrode';
    am = 'ActiveMaterial';
    co = 'Coating';
    bd = 'Binder';
    ca = 'ConductingAdditive';

    specVols = nan(3, 1);
    comps    = {am, bd, ca};
    eldes    = {ne, pe};

    for ielde = 1:numel(eldes)

        elde = eldes{ielde};

        for icomp = 1:numel(comps)

            comp = comps{icomp};
            rho  = paramobj.(elde).(co).(comp).density;
            mf   = paramobj.(elde).(co).(comp).massFraction;

            specVols(icomp) = mf/rho;
            assert(isfinite(specVols(icomp)));

        end

        paramobj.(elde).(co).volumeFractions = specVols ./ sum(specVols);
disp(elde)
        disp(        specVols ./ sum(specVols))
    end

    keyboard;

end
