classdef HighRateCalibration

    % Class to help perform high-rate calibration

    properties

        stdParams
        customParams
        customParamsSpec
        tag

    end

    methods


        function HRC = HighRateCalibration(simulatorSetup, tag)

            if nargin < 2
                HRC.tag = 'no-elyte-params';
            else
                HRC.tag = tag;
            end

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            itf   = 'Interface';
            sd    = 'SolidDiffusion';
            am    = 'ActiveMaterial';
            elyte = 'Electrolyte';
            sep   = 'Separator';
            rbc   = 'regionBruggemanCoefficients';

            eldes = {ne, pe};

            HRC.stdParams = [];

            % Setup params with standard getfun/setfun

            for ielde = 1:numel(eldes)

                elde = eldes{ielde};

                HRC.stdParams = addParameter(HRC.stdParams, ...
                                             simulatorSetup, ...
                                             'name'     , sprintf('%s_vsa', elde), ...
                                             'belongsTo', 'model'                , ...
                                             'boxLims'  , [1e4, 1e8]             , ...
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
                                             'boxLims', [0.1, 10], ...
                                             'scaling', 'linear', ...
                                             'getfun', @(model, ~) getEldeBruggeman(model), ...
                                             'setfun', @(model, ~, v) setEldeBruggeman(model, v), ...
                                             'location', {[{ne, co, 'bruggemanCoefficient'}; ...
                                                           {pe, co, 'bruggemanCoefficient'}]}); % location for printing

            switch HRC.tag
              case 'no-elyte-params'
                % Do nothing

              case {'one-elyte-param', 'one-elyte-param-finer'}

                HRC.stdParams = addParameter(HRC.stdParams, ...
                                             simulatorSetup, ...
                                             'name'     , 'elyte_bruggman', ...
                                             'belongsTo', 'model'                           , ...
                                             'boxLims'  , [0.1, 10]                         , ...
                                             'scaling'  , 'linear'                          , ...
                                             'location' , {elyte, 'bruggemanCoefficient'});

              case {'two-elyte-params', 'three-elyte-params'}

                HRC.customParamsSpec{end+1} = struct('name', 'elyte_bruggman', ...
                                                     'boxLims', [0.1, 10], ...
                                                     'scaling', 'linear', ...
                                                     'getfun', @(model, ~) getElyteBruggeman(model, tag), ...
                                                     'setfun', @(model, ~, v) setElyteBruggeman(model, v, tag), ...
                                                     'location', {[{elyte, rbc, ne}; ...
                                                                   {elyte, rbc, pe}; ...
                                                                   {elyte, rbc, sep}]}); % location for print

              otherwise
                error('Unknown tag: %s', HRC.tag);
            end

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
                for i = 1:size(vals, 1) % let the vals decide
                    jsonstruct = setfield(jsonstruct, locs{i,:}, vals(i));
                end
            end

        end

    end

end


function v = getEldeBruggeman(model)

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


function model = setEldeBruggeman(model, vals)

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


function v = getElyteBruggeman(model, tag)

    elyte = 'Electrolyte';
    w = model.(elyte).regionBruggemanCoefficients;

    switch tag
      case 'two-elyte-params'
        v = [w.NegativeElectrode;
             w.PositiveElectrode];
      case 'three-elyte-params'
        v = [w.NegativeElectrode;
             w.PositiveElectrode;
             w.Separator];
    end

end


function model = setElyteBruggeman(model, vals, tag)

    assert(~model.use_thermal);

    elyte = 'Electrolyte';
    ne    = 'NegativeElectrode';
    pe    = 'PositiveElectrode';
    sep   = 'Separator';

    nc    = model.(elyte).G.getNumberOfCells();
    tags  = model.(elyte).regionTags;
    bvals = model.(elyte).regionBruggemanCoefficients;
    bvals.(ne) = vals(1);
    bvals.(pe) = vals(2);

    bg = zeros(nc, 1);
    bg = subsetPlus(bg, bvals.(ne), (tags == 1));
    bg = subsetPlus(bg, bvals.(pe), (tags == 2));

    switch tag
      case 'three-elyte-params'
        bvals.(sep) = vals(3);
        bg = subsetPlus(bg, bvals.(sep), (tags == 3));
    end

    model.(elyte).bruggemanCoefficient = bg;
    model.(elyte).regionBruggemanCoefficients = bvals;

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
