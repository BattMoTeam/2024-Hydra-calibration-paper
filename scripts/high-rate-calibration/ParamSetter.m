classdef ParamSetter

    properties

    end

    methods

        function paramsetter = ParamSetter(model)

        end

        function vals = setFromVector(paramsetter, X)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';

            vals.(ne).vsa = X(1);
            vals.(pe).vsa = X(2);
            vals.(ne).bg  = X(3);
            vals.(pe).bg  = X(4);
            vals.(ne).D   = X(5);
            vals.(pe).D   = X(6);
            vals.(ne).ebg = X(7);
            vals.(pe).ebg = X(8);

        end

        function X = setToVector(paramsetter, vals)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';

            X(1, 1) = vals.(ne).vsa;
            X(2, 1) = vals.(pe).vsa;
            X(3, 1) = vals.(ne).bg;
            X(4, 1) = vals.(pe).bg;
            X(5, 1) = vals.(ne).D;
            X(6, 1) = vals.(pe).D;
            X(7, 1) = vals.(ne).ebg;
            X(8, 1) = vals.(pe).ebg;

        end

        function model = setValues(paramsetter, model, X)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            elyte   = 'Electrolyte';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';
            sep     = 'Separator';
            thermal = 'ThermalModel';

            vals = paramsetter.setFromVector(X);

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                bg = vals.(elde).bg;

                kappa = model.(elde).(co).electronicConductivity;
                vf    = model.(elde).(co).volumeFraction;

                model.(elde).(co).bruggemanCoefficient            = bg;
                model.(elde).(co).effectiveElectronicConductivity = kappa*vf^bg;

                model.(elde).(co).(am).(sd).referenceDiffusionCoefficient = vals.(elde).D;

                if model.use_thermal

                    emodel = model.(elde).(co);
                    thermalConductivity = 0;
                    for icomp = 1 : numel(emodel.compnames)
                        compname = emodel.compnames{icomp};
                        if ~isempty(emodel.(compname).thermalConductivity)
                            thermalConductivity = subsetPlus(thermalConductivity, (emodel.volumeFractions(icomp)).^bg*emodel.(compname).thermalConductivity, 1);
                        end
                    end
                    model.(elde).(co).thermalConductivity = thermalConductivity;

                    model.(elde).(co).effectiveThermalConductivity = (vf.^bg).*thermalConductivity;

                end

                model.(elde).(co).(am).(sd).volumetricSurfaceArea  = vals.(elde).vsa;
                model.(elde).(co).(am).(itf).volumetricSurfaceArea = vals.(elde).vsa;

            end

            nc    = model.(elyte).G.getNumberOfCells();
            tags  = model.(elyte).regionTags;
            bvals = model.(elyte).regionBruggemanCoefficients;
            bvals.(ne) = vals.(ne).ebg;
            bvals.(pe) = vals.(pe).ebg;
            model.(elyte).regionBruggemanCoefficients = bvals;

            bg = zeros(nc, 1);
            bg = subsetPlus(bg, bvals.NegativeElectrode, (tags == 1));
            bg = subsetPlus(bg, bvals.PositiveElectrode, (tags == 2));
            bg = subsetPlus(bg, bvals.Separator        , (tags == 3));

            model.(elyte).bruggemanCoefficient = bg;

            if model.use_thermal

                vf = model.(elyte).volumeFraction;
                model.(elyte).effectiveThermalConductivity = vf .^ bg .* model.(elyte).thermalConductivity;

            end


            if model.use_thermal

                G = model.G;
                nc = G.getNumberOfCells();

                hcond = zeros(nc, 1); % effective heat conductivity

                for ind = 1 : numel(eldes)

                    elde = eldes{ind};

                    if model.include_current_collectors

                        % The effecive and intrinsic thermal parameters for the current collector are the same.
                        cc_map   = model.(elde).(cc).G.mappings.cellmap;
                        cc_hcond = model.(elde).(cc).effectiveThermalConductivity;

                        hcond = subsetPlus(hcond, cc_hcond, cc_map);

                    end

                    % Effective parameters from the Electrode Active Component region.
                    co_map   = model.(elde).(co).G.mappings.cellmap;
                    co_hcond = model.(elde).(co).effectiveThermalConductivity;

                    hcond = subsetPlus(hcond, co_hcond, co_map);

                end

                % Electrolyte

                elyte_map   = model.(elyte).G.mappings.cellmap;
                elyte_hcond = model.(elyte).effectiveThermalConductivity;

                hcond = subsetPlus(hcond, elyte_hcond, elyte_map);

                % Separator

                sep_map   = model.(sep).G.mappings.cellmap;
                sep_hcond = model.(sep).effectiveThermalConductivity;

                hcond = subsetPlus(hcond, sep_hcond, sep_map);

                % Assign values

                model.(thermal).effectiveThermalConductivity = hcond;

            end

        end

        function X = getValues(paramsetter, model)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';

            eldes = {ne, pe};

            bvals = model.(elyte).regionBruggemanCoefficients;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                vals.(elde).bg  = model.(elde).(co).bruggemanCoefficient;
                vals.(elde).vsa = model.(elde).(co).(am).(sd).volumetricSurfaceArea;
                vals.(elde).D   = model.(elde).(co).(am).(sd).referenceDiffusionCoefficient;
                vals.(elde).ebg = bvals.(elde);

            end

            X = paramsetter.setToVector(vals);

        end

        function print(paramsetter, X)

            vals = paramsetter.setFromVector(X);

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            elyte = 'Electrolyte';

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                fprintf('%s Volumetric surface area : %g\n', elde, vals.(elde).vsa);
                fprintf('%s Bruggeman coefficient : %g\n'  , elde, vals.(elde).bg);
                fprintf('%s Diffusion coefficient : %g\n' , elde, vals.(elde).D);
            end

            fprintf('Electrolyte Bruggeman coefficient (negative region) : %g\n' , vals.(ne).ebg);
            fprintf('Electrolyte Bruggeman coefficient (positive region) : %g\n' , vals.(pe).ebg);

        end

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
