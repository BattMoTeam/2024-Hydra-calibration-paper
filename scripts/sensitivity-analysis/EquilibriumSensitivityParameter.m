classdef EquilibriumSensitivityParameter < ModelParameter

    % ModelParameter that keeps coupled EQC model and initial-state values
    % consistent when parameters are seeded for automatic differentiation.

    properties

        electrode
        parameterKind
        referenceLithiumTransfer
        npRatio

    end

    methods

        function parameter = EquilibriumSensitivityParameter(setup, varargin)

            opt = struct( ...
                'name', '', ...
                'electrode', '', ...
                'parameterKind', '', ...
                'boxLims', [], ...
                'referenceLithiumTransfer', [], ...
                'npRatio', []);
            opt = merge_options(opt, varargin{:});

            getfun = @(model, ~) ...
                EquilibriumSensitivityParameter.getValue( ...
                    model, opt.electrode, opt.parameterKind);
            unusedSetfun = @(model, ~, ~) model;

            parameter@ModelParameter( ...
                setup, ...
                'name', opt.name, ...
                'belongsTo', 'model', ...
                'boxLims', opt.boxLims, ...
                'scaling', 'linear', ...
                'location', {''}, ...
                'getfun', getfun, ...
                'setfun', unusedSetfun);

            parameter.electrode = opt.electrode;
            parameter.parameterKind = opt.parameterKind;
            parameter.referenceLithiumTransfer = ...
                opt.referenceLithiumTransfer;
            parameter.npRatio = opt.npRatio;

        end


        function setup = setParameter(parameter, setup, value)

            model = setup.model;

            switch parameter.parameterKind
              case 'theta100'
                model = EquilibriumSensitivityParameter.setTheta100( ...
                    model, parameter.electrode, value);
              case 'totalAmount'
                model = EquilibriumSensitivityParameter.setTotalAmount( ...
                    model, parameter.electrode, value);
              otherwise
                error('Unsupported EQC parameter kind: %s', ...
                      parameter.parameterKind);
            end

            model = EquilibriumSensitivityParameter.updateDerivedTheta0( ...
                model, parameter.referenceLithiumTransfer, parameter.npRatio);

            setup.model = model;
            setup.initstate = model.setupInitialState();

        end

    end

    methods(Static, Access = private)

        function value = getValue(model, electrode, parameterKind)

            co = 'Coating';
            am = 'ActiveMaterial';
            itf = 'Interface';

            switch parameterKind
              case 'theta100'
                value = model.(electrode).(co).(am).(itf) ...
                    .guestStoichiometry100;
              case 'totalAmount'
                values = EquilibriumSensitivityParameter.getElectrodeValues( ...
                    model, electrode);
                value = values.totalAmount;
              otherwise
                error('Unsupported EQC parameter kind: %s', parameterKind);
            end

        end


        function model = setTheta100(model, electrode, value)

            co = 'Coating';
            am = 'ActiveMaterial';
            itf = 'Interface';
            sd = 'SolidDiffusion';

            model.(electrode).(co).(am).(itf).guestStoichiometry100 = value;
            model.(electrode).(co).(am).(sd).guestStoichiometry100 = value;

        end


        function model = setTotalAmount(model, electrode, value)

            assert(~model.use_thermal, ...
                   'EQC volume-fraction sensitivities do not support thermal models.');

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            co = 'Coating';
            am = 'ActiveMaterial';
            sd = 'SolidDiffusion';
            elyte = 'Electrolyte';

            coating = model.(electrode).(co);
            activeMaterialFraction = ...
                coating.volumeFractions(coating.compInds.(am));
            amountPerVolumeFraction = ...
                sum(coating.G.getVolumes()) .* activeMaterialFraction .* ...
                coating.(am).Interface.saturationConcentration;
            volumeFraction = value ./ amountPerVolumeFraction;

            coating.volumeFraction = volumeFraction;
            coating.(am).volumeFraction = ...
                volumeFraction .* activeMaterialFraction;
            coating.(am).(sd).volumeFraction = ...
                volumeFraction .* activeMaterialFraction;
            model.(electrode).(co) = coating;

            electrolyteCells = zeros(model.G.getNumberOfCells(), 1);
            electrolyteCells(model.(elyte).G.mappings.cellmap) = ...
                (1:model.(elyte).G.getNumberOfCells())';
            regionCells = electrolyteCells( ...
                model.(electrode).(co).G.mappings.cellmap);

            electrolyteVolumeFraction = model.(elyte).volumeFraction;
            electrolyteVolumeFraction = subsasgnAD( ...
                electrolyteVolumeFraction, regionCells, ...
                1 - volumeFraction);
            model.(elyte).volumeFraction = electrolyteVolumeFraction;

            assert(any(strcmp(electrode, {ne, pe})), ...
                   'Unknown electrode: %s', electrode);

        end


        function model = updateDerivedTheta0( ...
                model, referenceLithiumTransfer, npRatio)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            neValues = EquilibriumSensitivityParameter.getElectrodeValues(model, ne);
            peValues = EquilibriumSensitivityParameter.getElectrodeValues(model, pe);

            peTheta0 = peValues.theta100 + ...
                referenceLithiumTransfer ./ peValues.totalAmount;
            neTheta0 = neValues.theta100 - ...
                npRatio .* referenceLithiumTransfer ./ neValues.totalAmount;

            model = EquilibriumSensitivityParameter.setTheta0(model, pe, peTheta0);
            model = EquilibriumSensitivityParameter.setTheta0(model, ne, neTheta0);

        end


        function model = setTheta0(model, electrode, value)

            co = 'Coating';
            am = 'ActiveMaterial';
            itf = 'Interface';
            sd = 'SolidDiffusion';

            model.(electrode).(co).(am).(itf).guestStoichiometry0 = value;
            model.(electrode).(co).(am).(sd).guestStoichiometry0 = value;

        end


        function values = getElectrodeValues(model, electrode)

            co = 'Coating';
            am = 'ActiveMaterial';
            itf = 'Interface';

            coating = model.(electrode).(co);
            interface = coating.(am).(itf);
            activeMaterialFraction = ...
                coating.volumeFractions(coating.compInds.(am));

            values = struct( ...
                'theta100', interface.guestStoichiometry100, ...
                'totalAmount', sum(coating.G.getVolumes()) .* ...
                    coating.volumeFraction .* activeMaterialFraction .* ...
                    interface.saturationConcentration);

        end

    end

end
