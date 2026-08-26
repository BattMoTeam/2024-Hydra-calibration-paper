classdef EquilibriumCalibrationSensitivity

    % Parameter catalog for adjoint sensitivities of the four independent
    % equilibrium-calibration (EQC) variables.

    properties

        params
        shortnames
        boxLims
        referenceLithiumTransfer
        npRatio

    end

    methods

        function EQC = EquilibriumCalibrationSensitivity(simulatorSetup)

            model = simulatorSetup.model;
            [referenceLithiumTransfer, npRatio] = ...
                EquilibriumCalibrationSensitivity.getCalibrationInvariants(model);

            EQC.referenceLithiumTransfer = referenceLithiumTransfer;
            EQC.npRatio = npRatio;

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            specs = [ ...
                EquilibriumCalibrationSensitivity.makeSpec( ...
                    model, ne, 'ne_theta100', 'theta100'), ...
                EquilibriumCalibrationSensitivity.makeSpec( ...
                    model, ne, 'ne_total_amount', 'totalAmount'), ...
                EquilibriumCalibrationSensitivity.makeSpec( ...
                    model, pe, 'pe_theta100', 'theta100'), ...
                EquilibriumCalibrationSensitivity.makeSpec( ...
                    model, pe, 'pe_total_amount', 'totalAmount')];

            EQC.shortnames = reshape({specs.shortname}, [], 1);
            EQC.boxLims = vertcat(specs.boxLims);
            EQC.params = cell(numel(specs), 1);

            for index = 1:numel(specs)
                spec = specs(index);
                EQC.params{index} = EquilibriumSensitivityParameter( ...
                    simulatorSetup, ...
                    'name', spec.shortname, ...
                    'electrode', spec.electrode, ...
                    'parameterKind', spec.parameterKind, ...
                    'boxLims', spec.boxLims, ...
                    'referenceLithiumTransfer', referenceLithiumTransfer, ...
                    'npRatio', npRatio);
            end

        end


        function params = getParams(EQC)

            params = EQC.params;

        end

    end

    methods(Static, Access = private)

        function spec = makeSpec(model, electrode, shortname, parameterKind)

            switch parameterKind
              case 'theta100'
                boxLims = [0, 1];
              case 'totalAmount'
                values = EquilibriumCalibrationSensitivity.getElectrodeValues( ...
                    model, electrode);
                boxLims = values.totalAmount .* [0.1, 10];
              otherwise
                error('Unsupported EQC parameter kind: %s', parameterKind);
            end

            spec = struct( ...
                'shortname', shortname, ...
                'electrode', electrode, ...
                'parameterKind', parameterKind, ...
                'boxLims', boxLims);

        end


        function [referenceLithiumTransfer, npRatio] = ...
                getCalibrationInvariants(model)

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';

            neValues = EquilibriumCalibrationSensitivity.getElectrodeValues(model, ne);
            peValues = EquilibriumCalibrationSensitivity.getElectrodeValues(model, pe);

            peWindow = peValues.theta0 - peValues.theta100;
            neWindow = neValues.theta100 - neValues.theta0;

            assert(peWindow > 0 && neWindow > 0, ...
                   'Expected positive calibrated stoichiometric windows.');

            referenceLithiumTransfer = peValues.totalAmount .* peWindow;
            npRatio = neValues.totalAmount .* neWindow ./ ...
                      referenceLithiumTransfer;

            assert(all([referenceLithiumTransfer, npRatio] > 0), ...
                   'The EQC calibration invariants must be positive.');

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
                'theta0', interface.guestStoichiometry0, ...
                'totalAmount', sum(coating.G.getVolumes()) .* ...
                    coating.volumeFraction .* activeMaterialFraction .* ...
                    interface.saturationConcentration);

        end

    end

end
