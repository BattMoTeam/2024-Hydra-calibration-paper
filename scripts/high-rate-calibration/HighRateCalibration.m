classdef HighRateCalibration

    % Class to help perform high-rate calibration

    properties

        params
        parameterSpecs
        shortnames
        boxLims
        useRegionBruggemanCoefficients

    end

    methods

        function HRC = HighRateCalibration(simulatorSetup, varargin)

            opt = struct('shortnames', []);
            opt = merge_options(opt, varargin{:});

            elyte = 'Electrolyte';
            HRC.useRegionBruggemanCoefficients = ...
                simulatorSetup.model.(elyte).useRegionBruggemanCoefficients;

            specs = HRC.parameterCatalog();
            availableShortnames = {specs.shortname};

            if ~isempty(opt.shortnames)
                selected = HighRateCalibration.normalizeShortnames(opt.shortnames);

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

            HRC.parameterSpecs = reshape(specs, [], 1);
            HRC.shortnames = reshape({specs.shortname}, [], 1);
            HRC.boxLims = vertcat(specs.boxLim);
            HRC.params = cell(numel(specs), 1);

            for index = 1:numel(specs)
                spec = specs(index);

                if isempty(spec.getfun)
                    HRC.params{index} = ModelParameter( ...
                        simulatorSetup, ...
                        'name'     , spec.shortname, ...
                        'belongsTo', 'model', ...
                        'boxLims'  , spec.boxLim, ...
                        'scaling'  , spec.scaling, ...
                        'location' , spec.location);
                else
                    HRC.params{index} = ModelParameter( ...
                        simulatorSetup, ...
                        'name'     , spec.shortname, ...
                        'belongsTo', 'model', ...
                        'boxLims'  , spec.boxLim, ...
                        'scaling'  , spec.scaling, ...
                        'location' , {''}, ...
                        'getfun'   , spec.getfun, ...
                        'setfun'   , spec.setfun);
                end
            end

        end


        function specs = parameterCatalog(HRC)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            itf   = 'Interface';
            sd    = 'SolidDiffusion';
            am    = 'ActiveMaterial';
            elyte = 'Electrolyte';
            sep   = 'Separator';
            rbc   = 'regionBruggemanCoefficients';

            specs = struct('shortname', {}, ...
                           'location' , {}, ...
                           'boxLim'   , {}, ...
                           'scaling'  , {}, ...
                           'getfun'   , {}, ...
                           'setfun'   , {});

            specs(end + 1) = struct( ...
                'shortname', 'ne_vsa', ...
                'location' , {{ne, co, am, itf, 'volumetricSurfaceArea'}}, ...
                'boxLim'   , [1e3, 1e8], ...
                'scaling'  , 'log', ...
                'getfun'   , @(model, ~) getVsa(model, ne), ...
                'setfun'   , @(model, ~, v) setVsa(model, v, ne));

            specs(end + 1) = struct( ...
                'shortname', 'pe_vsa', ...
                'location' , {{pe, co, am, itf, 'volumetricSurfaceArea'}}, ...
                'boxLim'   , [1e4, 1e8], ...
                'scaling'  , 'log', ...
                'getfun'   , @(model, ~) getVsa(model, pe), ...
                'setfun'   , @(model, ~, v) setVsa(model, v, pe));

            specs(end + 1) = struct( ...
                'shortname', 'ne_D', ...
                'location' , {{ne, co, am, sd, 'referenceDiffusionCoefficient'}}, ...
                'boxLim'   , [1e-15, 1e-11], ...
                'scaling'  , 'log', ...
                'getfun'   , [], ...
                'setfun'   , []);

            specs(end + 1) = struct( ...
                'shortname', 'pe_D', ...
                'location' , {{pe, co, am, sd, 'referenceDiffusionCoefficient'}}, ...
                'boxLim'   , [1e-15, 1e-11], ...
                'scaling'  , 'log', ...
                'getfun'   , [], ...
                'setfun'   , []);

            specs(end + 1) = struct( ...
                'shortname', 'ne_bg', ...
                'location' , {{ne, co, 'bruggemanCoefficient'}}, ...
                'boxLim'   , [0.1, 10], ...
                'scaling'  , 'linear', ...
                'getfun'   , @(model, ~) getEldeBruggeman(model, ne), ...
                'setfun'   , @(model, ~, v) setEldeBruggeman(model, v, ne));

            specs(end + 1) = struct( ...
                'shortname', 'pe_bg', ...
                'location' , {{pe, co, 'bruggemanCoefficient'}}, ...
                'boxLim'   , [0.1, 10], ...
                'scaling'  , 'linear', ...
                'getfun'   , @(model, ~) getEldeBruggeman(model, pe), ...
                'setfun'   , @(model, ~, v) setEldeBruggeman(model, v, pe));

            if HRC.useRegionBruggemanCoefficients
                specs(end + 1) = struct( ...
                    'shortname', 'elyte_bg_ne', ...
                    'location' , {{elyte, rbc, ne}}, ...
                    'boxLim'   , [0.1, 10], ...
                    'scaling'  , 'linear', ...
                    'getfun'   , @(model, ~) getElyteRegionBruggeman(model, ne), ...
                    'setfun'   , @(model, ~, v) setElyteRegionBruggeman(model, v, ne));

                specs(end + 1) = struct( ...
                    'shortname', 'elyte_bg_pe', ...
                    'location' , {{elyte, rbc, pe}}, ...
                    'boxLim'   , [0.1, 10], ...
                    'scaling'  , 'linear', ...
                    'getfun'   , @(model, ~) getElyteRegionBruggeman(model, pe), ...
                    'setfun'   , @(model, ~, v) setElyteRegionBruggeman(model, v, pe));

                specs(end + 1) = struct( ...
                    'shortname', 'elyte_bg_sep', ...
                    'location' , {{elyte, rbc, sep}}, ...
                    'boxLim'   , [0.1, 10], ...
                    'scaling'  , 'linear', ...
                    'getfun'   , @(model, ~) getElyteRegionBruggeman(model, sep), ...
                    'setfun'   , @(model, ~, v) setElyteRegionBruggeman(model, v, sep));

                specs(end + 1) = struct( ...
                    'shortname', 'elyte_bgfactor', ...
                    'location' , {{elyte, 'bgfactor'}}, ...
                    'boxLim'   , [0.1, 10], ...
                    'scaling'  , 'linear', ...
                    'getfun'   , @(model, ~) getElyteBgfactor(model), ...
                    'setfun'   , @(model, ~, v) setElyteBgfactor(model, v));
            else
                specs(end + 1) = struct( ...
                    'shortname', 'elyte_bg', ...
                    'location' , {{elyte, 'bruggemanCoefficient'}}, ...
                    'boxLim'   , [0.1, 10], ...
                    'scaling'  , 'linear', ...
                    'getfun'   , [], ...
                    'setfun'   , []);
            end

        end


        function specs = selectedParameterCatalog(HRC)

            specs = HRC.parameterSpecs;

        end


        function locs = locations(HRC)

            specs = HRC.selectedParameterCatalog();
            locs = reshape({specs.location}, [], 1);

        end


        function params = getParams(HRC)

            params = HRC.params;

        end


        function jsonstruct = export(HRC, setup)

            jsonstruct = struct();

            for index = 1:numel(HRC.params)
                loc = HRC.parameterSpecs(index).location;
                value = HRC.params{index}.getParameterValue(setup);
                assert(isscalar(value), 'Expected scalar value for %s.', HRC.shortnames{index});
                jsonstruct = setfield(jsonstruct, loc{:}, value);
            end

        end


        function printBoxLims(HRC)

            disp('boxLims:');
            tbl = table(HRC.shortnames, HRC.boxLims(:, 1), HRC.boxLims(:, 2), ...
                        'VariableNames', {'Shortname', 'LowerLimit', 'UpperLimit'});
            disp(tbl);

        end

    end

    methods (Static)

        function shortnames = normalizeShortnames(shortnames)

            if ischar(shortnames)
                shortnames = {shortnames};
            elseif isstring(shortnames)
                shortnames = cellstr(shortnames(:));
            elseif iscell(shortnames)
                shortnames = shortnames(:);
                shortnames = cellfun(@char, shortnames, 'UniformOutput', false);
            else
                error('''shortnames'' must be a char, string array, or cell array of char.');
            end

        end

    end

end


function value = getVsa(model, electrode)

    co  = 'Coating';
    am  = 'ActiveMaterial';
    itf = 'Interface';

    value = model.(electrode).(co).(am).(itf).volumetricSurfaceArea;

end


function model = setVsa(model, value, electrode)

    co  = 'Coating';
    am  = 'ActiveMaterial';
    itf = 'Interface';
    sd  = 'SolidDiffusion';

    model.(electrode).(co).(am).(itf).volumetricSurfaceArea = value;
    model.(electrode).(co).(am).(sd).volumetricSurfaceArea = value;

end


function value = getEldeBruggeman(model, electrode)

    co = 'Coating';
    value = model.(electrode).(co).bruggemanCoefficient;

end


function model = setEldeBruggeman(model, value, electrode)

    assert(~model.use_thermal);

    co = 'Coating';
    model.(electrode).(co).bruggemanCoefficient = value;

    kappa = model.(electrode).(co).electronicConductivity;
    volumeFraction = model.(electrode).(co).volumeFraction;
    model.(electrode).(co).effectiveElectronicConductivity = ...
        kappa .* volumeFraction .^ value;

end


function value = getElyteRegionBruggeman(model, region)

    elyte = 'Electrolyte';
    value = model.(elyte).regionBruggemanCoefficients.(region);

end


function model = setElyteRegionBruggeman(model, value, region)

    assert(~model.use_thermal);

    elyte = 'Electrolyte';
    model.(elyte).regionBruggemanCoefficients.(region) = value;
    model = updateElyteBruggemanCoefficient(model);

end


function value = getElyteBgfactor(model)

    elyte = 'Electrolyte';
    value = model.(elyte).bgfactor;

end


function model = setElyteBgfactor(model, value)

    assert(~model.use_thermal);

    model.Electrolyte.bgfactor = value;
    % model = updateElyteBruggemanCoefficient(model);

end


function model = updateElyteBruggemanCoefficient(model)

    elyte = 'Electrolyte';
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    sep   = 'Separator';

    numberOfCells = model.(elyte).G.getNumberOfCells();
    tags = model.(elyte).regionTags;
    regionValues = model.(elyte).regionBruggemanCoefficients;
    % bgfactor = model.(elyte).bgfactor;

    % bg = zeros(numberOfCells, 1);
    % bg = subsetPlus(bg, bgfactor .* regionValues.(ne), tags == 1);
    % bg = subsetPlus(bg, bgfactor .* regionValues.(pe), tags == 2);
    % bg = subsetPlus(bg, bgfactor .* regionValues.(sep), tags == 3);

    bg = zeros(numberOfCells, 1);
    bg = subsetPlus(bg, regionValues.(ne), tags == 1);
    bg = subsetPlus(bg, regionValues.(pe), tags == 2);
    bg = subsetPlus(bg, regionValues.(sep), tags == 3);

    model.(elyte).bruggemanCoefficient = bg;

end


%{
  Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
  and SINTEF Digital, Mathematics & Cybernetics.

  This file is part of The Battery Modeling Toolbox BattMo

  BattMo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BattMo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
