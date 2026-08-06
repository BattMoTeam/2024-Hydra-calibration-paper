classdef ParamSetter

    properties
        include_current_collectors
        use_equivalent_eff_cond
        useRegionBruggemanCoefficients

        shortnames
        boxLims
    end

    methods

        function paramsetter = ParamSetter(varargin)

            opt = struct('useRegionBruggemanCoefficients'  , false, ...
                         'shortnames'                      , []   );
            opt = merge_options(opt, varargin{:});

            paramsetter.useRegionBruggemanCoefficients = opt.useRegionBruggemanCoefficients;

            specs = paramsetter.parameterCatalog();
            availableShortnames = {specs.shortname};

            if isempty(opt.shortnames)
                selected = availableShortnames;
            else
                selected = ParamSetter.normalizeShortnames(opt.shortnames);

                if numel(unique(selected)) ~= numel(selected)
                    error('The ''shortnames'' list contains duplicates.');
                end

                [tf, idx] = ismember(selected, availableShortnames);
                if ~all(tf)
                    bad = selected(~tf);
                    error('Unknown or unavailable shortnames: %s', strjoin(bad, ', '));
                end

                specs = specs(idx); % preserve user-specified order
            end

            assert(~isempty(specs), 'At least one parameter must be selected.');

            paramsetter.shortnames = reshape({specs.shortname}, [], 1);
            paramsetter.boxLims = vertcat(specs.boxLim);

        end

        function specs = parameterCatalog(paramsetter)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            elyte = 'Electrolyte';
            sep   = 'Separator';

            specs = struct('shortname', {}, 'location', {}, 'boxLim', {});

            specs(end+1) = struct('shortname', 'SEI_L0', ...
                                  'location' , {{ne, co, am, itf, 'SEIlengthRef'}}, ...
                                  'boxLim'   , [1e-10, 1e-6]);

            specs(end+1) = struct('shortname', 'SEI_D', ...
                                  'location' , {{ne, co, am, itf, 'SEIelectronicDiffusionCoefficient'}}, ...
                                  'boxLim'   , [1e-16, 1e-8]);

            specs(end+1) = struct('shortname', 'SEI_kappa', ...
                                  'location' , {{ne, co, am, itf, 'SEIionicConductivity'}}, ...
                                  'boxLim'   , [1e-9, 1e-1]);

            specs(end+1) = struct('shortname', 'ne_j0', ...
                                  'location' , {{ne, co, am, itf, 'reactionRateConstant'}}, ...
                                  'boxLim'   , [1e-13, 1e-6]);

            specs(end+1) = struct('shortname', 'pe_j0', ...
                                  'location' , {{pe, co, am, itf, 'reactionRateConstant'}}, ...
                                  'boxLim'   , [1e-13, 1e-9]);

            specs(end+1) = struct('shortname', 'ne_vsa', ...
                                  'location' , {{ne, co, am, itf, 'volumetricSurfaceArea'}}, ...
                                  'boxLim'   , [1e4, 1e9]);

            specs(end+1) = struct('shortname', 'pe_vsa', ...
                                  'location' , {{pe, co, am, itf, 'volumetricSurfaceArea'}}, ...
                                  'boxLim'   , [1e5, 1e11]);

            specs(end+1) = struct('shortname', 'ne_bg', ...
                                  'location' , {{ne, co, 'bruggemanCoefficient'}}, ...
                                  'boxLim'   , [1e-10, 20]);

            specs(end+1) = struct('shortname', 'pe_bg', ...
                                  'location' , {{pe, co, 'bruggemanCoefficient'}}, ...
                                  'boxLim'   , [1e-10, 20]);

            specs(end+1) = struct('shortname', 'ne_D', ...
                                  'location' , {{ne, co, am, sd, 'referenceDiffusionCoefficient'}}, ...
                                  'boxLim'   , [1e-15, 1e-9]);

            specs(end+1) = struct('shortname', 'pe_D', ...
                                  'location' , {{pe, co, am, sd, 'referenceDiffusionCoefficient'}}, ...
                                  'boxLim'   , [1e-15, 1e-8]);

            if paramsetter.useRegionBruggemanCoefficients
                specs(end+1) = struct('shortname', 'elyte_bg_ne', ...
                                      'location' , {{elyte, 'regionBruggemanCoefficients', ne}}, ...
                                      'boxLim'   , [1e-10, 10.0]);

                specs(end+1) = struct('shortname', 'elyte_bg_pe', ...
                                      'location' , {{elyte, 'regionBruggemanCoefficients', pe}}, ...
                                      'boxLim'   , [1e-10, 10.0]);

                specs(end+1) = struct('shortname', 'elyte_bg_sep', ...
                                      'location' , {{elyte, 'regionBruggemanCoefficients', sep}}, ...
                                      'boxLim'   , [1e-10, 10.0]);
            else
                specs(end+1) = struct('shortname', 'elyte_bg', ...
                                      'location' , {{elyte, 'bruggemanCoefficient'}}, ...
                                      'boxLim'   , [1e-10, 10.0]);
            end

        end

        function specs = selectedParameterCatalog(paramsetter)
            allSpecs = paramsetter.parameterCatalog();
            allShortnames = {allSpecs.shortname};

            [tf, idx] = ismember(paramsetter.shortnames, allShortnames);
            assert(all(tf), 'Selected shortnames are inconsistent with the current options.');

            specs = allSpecs(idx);
        end

        function locs = locations(paramsetter)
            specs = paramsetter.selectedParameterCatalog();
            locs = reshape({specs.location}, [], 1);

            assert(numel(locs) == size(paramsetter.boxLims,1), ...
                   'Number of locations and box limits must match');
        end

        function isValid = validate(paramsetter, model)

            locs = paramsetter.locations();
            isValid = false(numel(locs), 1);

            for k = 1:numel(locs)
                loc = locs{k};
                subloc = loc(1:end-1);
                submodel = getfield(model, subloc{:});
                if isprop(submodel, loc{end}) || isfield(submodel, loc{end})
                    isValid(k) = true;
                end
            end

            if ~all(isValid)
                sn = paramsetter.shortnames();
                tbl = table(sn, isValid);
                disp(tbl);
            end

            assert(all(isValid), 'The model does not contain all the required fields');

        end

        function s = shortlocs(paramsetter)

            sn = paramsetter.shortnames;
            s = cell(numel(sn), 1);

            for k = 1:numel(sn)
                s{k} = regexp(sn{k}, '_', 'split');
            end

        end

        function paramsetter = setupRelBox(paramsetter, relsize, X)

            if isscalar(relsize)
                relsize = relsize*ones(size(X));
            end

            boxLims = paramsetter.boxLims;
            for k = 1:numel(X)
                boxLims(k, 1) = X(k)*(1 - relsize(k));
                boxLims(k, 2) = X(k)*(1 + relsize(k));
            end

            paramsetter.boxLims = boxLims;

        end

        function vals = setFromVector(paramsetter, X)

            slocs = paramsetter.shortlocs();
            vals = struct();

            for k = 1:numel(slocs)
                s = slocs{k};
                vals = setfield(vals, s{:}, X(k));
            end

        end

        function X = setToVector(paramsetter, vals)

            locs = paramsetter.shortlocs();
            X = nan(numel(locs), 1);

            for k = 1:numel(locs)
                loc = locs{k};
                X(k) = getfield(vals, loc{:});
            end

        end

        function model = setValues(paramsetter, model, X)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';

            % Automatic update of explicitly selected parameters
            locs = paramsetter.locations();
            for k = 1:numel(locs)
                loc = locs{k};
                model = setfield(model, loc{:}, X(k));
            end

            % Selected parameter names
            sn = paramsetter.shortnames;

            if any(strcmp(sn, 'ne_bg'))
                bg    = model.(ne).(co).bruggemanCoefficient;
                kappa = model.(ne).(co).electronicConductivity;
                vf    = model.(ne).(co).volumeFraction;
                model.(ne).(co).effectiveElectronicConductivity = kappa*vf^bg;
            end

            if any(strcmp(sn, 'pe_bg'))
                bg    = model.(pe).(co).bruggemanCoefficient;
                kappa = model.(pe).(co).electronicConductivity;
                vf    = model.(pe).(co).volumeFraction;
                model.(pe).(co).effectiveElectronicConductivity = kappa*vf^bg;
            end

            if any(strcmp(sn, 'ne_vsa'))
                model.(ne).(co).(am).(sd).volumetricSurfaceArea = ...
                    model.(ne).(co).(am).(itf).volumetricSurfaceArea;
            end

            if any(strcmp(sn, 'pe_vsa'))
                model.(pe).(co).(am).(sd).volumetricSurfaceArea = ...
                    model.(pe).(co).(am).(itf).volumetricSurfaceArea;
            end

            if model.(elyte).useRegionBruggemanCoefficients

                elyteBgShortnames = {'elyte_bg_ne', 'elyte_bg_pe', 'elyte_bg_sep'};

                if any(ismember(sn, elyteBgShortnames))

                    nc = model.(elyte).G.getNumberOfCells();
                    bg = zeros(nc, 1);

                    tagmap = struct('NegativeElectrode', 1, ...
                                    'PositiveElectrode', 2, ...
                                    'Separator', 3);
                    tags = model.(elyte).regionTags;

                    regionNames = {'NegativeElectrode', 'PositiveElectrode', 'Separator'};

                    for i = 1:numel(regionNames)
                        region = regionNames{i};
                        bval = model.(elyte).regionBruggemanCoefficients.(region);
                        bg = subsetPlus(bg, bval, (tags == tagmap.(region)));
                    end

                    model.(elyte).bruggemanCoefficient = bg;
                end

            end

        end

        function X = getValues(paramsetter, model)

            locs  = paramsetter.locations();
            short = paramsetter.shortlocs();

            vals = struct();

            for k = 1:numel(locs)
                loc = locs{k};
                mval = getfield(model, loc{:});
                s = short{k};
                vals = setfield(vals, s{:}, mval);
            end

            X = paramsetter.setToVector(vals);

        end

        function print(paramsetter, X)

            printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));
            vals = paramsetter.setFromVector(X);
            printer(vals);

        end

        function printBoxLims(paramsetter)

            tbl = table(paramsetter.shortnames, paramsetter.boxLims(:,1), paramsetter.boxLims(:,2), ...
                        'VariableNames', {'Shortname', 'LowerLimit', 'UpperLimit'});
            disp(tbl);

        end

    end

    methods (Static)

        function sn = normalizeShortnames(shortnames)

            if ischar(shortnames)
                sn = {shortnames};
            elseif isstring(shortnames)
                sn = cellstr(shortnames(:));
            elseif iscell(shortnames)
                sn = shortnames(:);
                sn = cellfun(@char, sn, 'UniformOutput', false);
            else
                error('''shortnames'' must be a char, string array, or cell array of char.');
            end

        end

    end

end
