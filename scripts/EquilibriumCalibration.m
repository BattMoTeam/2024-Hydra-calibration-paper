classdef EquilibriumCalibration

% Class to help perform the low-rate calibration under equilibrium
% assumptions

% Parameters:
% Y(1): guestStoichiometry100 cathode
% Y(2): alpha cathode (alpha = V*volumeFraction*cmax)
% Y(3): guestStoichiometry100 anode
% Y(4): alpha anode (alpha = V*volumeFraction*cmax)

    properties

        F
        time
        E
        totalTime
        I
        DRate
        model

    end

    methods

        function EC = EquilibriumCalibration(model, expdata)

            EC.F         = PhysicalConstants().F;
            EC.time      = expdata.time(:);
            EC.E         = expdata.E;
            EC.I         = expdata.I;
            EC.DRate     = expdata.DRate;
            EC.totalTime = EC.time(end);
            EC.model     = model;

            assert(isscalar(EC.I));

        end


        function Y = getDefaultValue(EC)

            pe = 'PositiveElectrode';
            ne = 'NegativeElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            co  = 'Coating';

            Y = nan(4, 1);

            Y(1) = EC.model.(pe).(co).(am).(itf).guestStoichiometry100;
            Y(2) = computeAlpha(EC.model.(pe).(co), EC.model.(pe).(co).compInds.(am));

            Y(3) = EC.model.(ne).(co).(am).(itf).guestStoichiometry100;
            Y(4) = computeAlpha(EC.model.(ne).(co), EC.model.(ne).(co).compInds.(am));

        end


        function [z, dz] = objective(EC, X)

            [fexp, fcomp] = EC.setupFunctions();
            X = initVariablesADI(X);

            fdiff = (fcomp(EC.time, X) - fexp(EC.time)).^2;

            g = 0.5 * sum(diff(EC.time) .* (fdiff(1:end-1) + fdiff(2:end)));
            z = g.val;
            dz = g.jac{1}';

        end


        function [fexp, fcomp] = setupFunctions(EC)

            fexp = @(t) EC.experimentalF(t);
            fcomp = @(t, X) EC.computeF(t, X);

        end


        function c = conc(EC, t, sgn, theta, alpha)

            c = theta + t * ((sgn * EC.I) ./ (EC.F * alpha));

        end


        function f = computeF(EC, t, X)

            pe  = 'PositiveElectrode';
            ne  = 'NegativeElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            co  = 'Coating';

            T = [];

            theta = X(1);
            alpha = X(2);
            c = EC.conc(t, 1, theta, alpha);
            f = EC.model.(pe).(co).(am).(itf).computeOCPFunc(c, T, 1);

            theta = X(3);
            alpha = X(4);
            c = EC.conc(t, -1, theta, alpha);
            f = f - EC.model.(ne).(co).(am).(itf).computeOCPFunc(c, T, 1);

        end


        function f = experimentalF(EC, t)

            f = interp1(EC.time, EC.E, t);

        end


        function jsonstruct = extractAlpha(EC, paramobj, Yopt)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            co  = 'Coating';

            jsonstruct = struct();
            Xopt       = toStruct(Yopt);
            eldes      = {ne, pe};

            for ielde = 1:numel(eldes)
                elde = eldes{ielde};

                [theta0, theta100] = EC.calculateGuestStoichiometry(Xopt, elde);
                jsonstruct.(elde).(co).(am).(itf).guestStoichiometry0   = theta0;
                jsonstruct.(elde).(co).(am).(itf).guestStoichiometry100 = theta100;

                indam = EC.model.(elde).(co).compInds.(am);
                vf    = calculateVF(Xopt, paramobj, elde, indam);
                jsonstruct.(elde).(co).volumeFraction = vf;

            end

        end


        function [theta0, theta100] = calculateGuestStoichiometry(EC, Xopt, elde)

            theta100 = Xopt.(elde).theta100;

            if theta100 < 0 || theta100 > 1
                error('guestStoichiometry100 for %s is not between 0 and 1', elde);
            end

            switch elde
              case 'NegativeElectrode'
                sgn = -1;
              case 'PositiveElectrode'
                sgn = 1;
              otherwise
                error('Unknown electrode %g', elde);
            end

            theta0 = EC.conc(EC.totalTime, sgn, theta100, Xopt.(elde).alpha);

            if theta0 < 0 || theta0 > 1
                error('guestStoichiometry0 for %s is not between 0 and 1', elde);
            end

        end



    end


end


function Y = toStruct(X)

    pe = 'PositiveElectrode';
    ne = 'NegativeElectrode';

    Y.(pe) = struct('theta100', X(1), ...
                    'alpha'   , X(2));
    Y.(ne) = struct('theta100', X(3), ...
                    'alpha'   , X(4));

end


function alpha = computeAlpha(coating, indam)

    am  = 'ActiveMaterial';
    itf = 'Interface';

    vol  = sum(coating.G.getVolumes());
    vf   = coating.volumeFraction;
    vfam = coating.volumeFractions(indam);
    cmax = coating.(am).(itf).saturationConcentration;

    alpha = vol * vf * vfam * cmax;

end


function pf = getPropFactor(Xopt, paramobj, elde, indam)

    co = 'Coating';

    alphaOpt = Xopt.(elde).alpha;
    alpha    = computeAlpha(paramobj.(elde).(co), indam);
    pf       = alphaOpt / alpha;

end


function vf = calculateVF(Xopt, paramobj, elde, indam)

    co = 'Coating';

    p  = getPropFactor(Xopt, paramobj, elde, indam);
    vf = p * paramobj.(elde).(co).volumeFraction;

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
