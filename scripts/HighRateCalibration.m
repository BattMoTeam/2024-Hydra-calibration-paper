classdef HighRateCalibration

% Class to help perform high-rate calibration

    properties

        stdParams
        customParams
        customParamsSpec

    end

    methods


        function HRC = HighRateCalibration(simulatorSetup)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
            am  = 'ActiveMaterial';

            eldes = {ne, pe};

            HRC.stdParams = [];

            % Setup params with standard getfun/setfun

            for ielde = 1:numel(eldes)

                elde = eldes{ielde};

                HRC.stdParams = addParameter(HRC.stdParams, ...
                                             simulatorSetup, ...
                                             'name'     , sprintf('%s_vsa', elde), ...
                                             'belongsTo', 'model'                , ...
                                             'boxLims'  , [1e4, 1e7]             , ...
                                             'location' , {elde, co, am, itf, 'volumetricSurfaceArea'});

                HRC.stdParams = addParameter(HRC.stdParams, ...
                                             simulatorSetup, ...
                                             'name'     , sprintf('%s_D0', elde), ...
                                             'belongsTo', 'model'               , ...
                                             'boxLims'  , [1e-15, 1e-11]        , ...
                                             'scaling'  , 'log'                 , ...
                                             'location' , {elde, co, am, sd, 'referenceDiffusionCoefficient'});
            end

            % Setup params with custom getfun/setfun
            HRC.customParamsSpec{1} = struct('name', 'eldes_bruggeman', ...
                                             'boxLims', [1.5, 20], ...
                                             'scaling', 'linear', ...
                                             'getfun', @(model, ~) getBruggeman(model), ...
                                             'setfun', @(model, ~, v) setBruggeman(model, v), ...
                                             'location', {[{ne, co, 'bruggemanCoefficient'}; ...
                                                           {pe, co, 'bruggemanCoefficient'}]});

            % Convert spec to ModelParameter instances
            HRC.customParams = cell(numel(HRC.customParamsSpec), 1);

            for k = 1:numel(HRC.customParamsSpec)

                spec = HRC.customParamsSpec{k};

                HRC.customParams{k} = ModelParameter(simulatorSetup, ...
                                                     'name', spec.name, ...
                                                     'belongsTo', 'model', ...
                                                     'boxLims', spec.boxLims, ...
                                                     'scaling', spec.scaling, ...
                                                     'location', {''}, ...
                                                     'getfun', spec.getfun, ...
                                                     'setfun', spec.setfun);
            end

            HRC.stdParams = reshape(HRC.stdParams, [], 1);
            HRC.customParamsSpec = reshape(HRC.customParamsSpec, [], 1);
            HRC.customParams = reshape(HRC.customParams, [], 1);

        end


        function params = getParams(HRC)

            params = [HRC.stdParams; HRC.customParams];

        end


        function jsonstruct = export(HRC, setup)

            % Standard params
            locs_std = cellfun(@(p) p.location, HRC.stdParams, 'uniformoutput', false);
            vals_std = cellfun(@(p) p.getParameterValue(setup), HRC.stdParams);

            % Custom params
            locs_custom = cellfun(@(p) p.location, HRC.customParamsSpec, 'uniformoutput', false);
            vals_custom = cellfun(@(p) p.getParameterValue(setup), HRC.customParams, 'uniformoutput', false);

            jsonstruct = struct();

            for k = 1:numel(locs_std)
                loc = locs_std{k};
                jsonstruct = setfield(jsonstruct, loc{:}, vals_std(k));
            end

            for k = 1:numel(locs_custom)
                locs = locs_custom{k};
                vals = vals_custom{k};
                for i = 1:size(locs, 1)
                    jsonstruct = setfield(jsonstruct, locs{i,:}, vals(i));
                end
            end

        end

    end

end


function v = getBruggeman(model)

    ne = 'NegativeElectrode';
    pe = 'PositiveElectrode';
    co = 'Coating';
    eldes = {ne, pe};

    v = nan(numel(eldes), 1);

    for ielde = 1:numel(eldes)
        elde = eldes{ielde};
        v(ielde) = model.(elde).(co).bruggemanCoefficient;
    end

end


function model = setBruggeman(model, vals)

    assert(~model.use_thermal);

    ne = 'NegativeElectrode';
    pe = 'PositiveElectrode';
    co = 'Coating';
    eldes = {ne, pe};

    for ielde = 1:numel(eldes)
        elde = eldes{ielde};

        % Set value
        bg = vals(ielde);
        model.(elde).(co).bruggemanCoefficient = bg;

        % Set dependencies
        kappa = model.(elde).(co).electronicConductivity;
        vf = model.(elde).(co).volumeFraction;
        model.(elde).(co).effectiveElectronicConductivity = kappa*vf^bg;

    end

end


%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
