function D = computeDcathodeH0b(soc)

    lnmo = [
        0.14	1.3026E-14
        0.185263157894737	1.34825171070304E-14
        0.230526315789474	1.5093434208444E-14
        0.275789473684211	1.4773079647521E-14
        0.321052631578947	1.44616508830575E-14
        0.366315789473684	6.79590491166938E-15
        0.411578947368421	4.53414466608471E-15
        0.456842105263158	4.43405350345279E-15
        0.502105263157895	4.33396234082088E-15
        0.547368421052632	4.23387117818896E-15
        0.592631578947368	4.13378001555705E-15
        0.637894736842105	4.03368885292513E-15
        0.683157894736842	4.03302884223173E-15
        0.728421052631579	4.10328663788827E-15
        0.773684210526316	4.17354443354481E-15
        0.818947368421053	4.24380222920136E-15
        0.864210526315789	4.3140600248579E-15
        0.909473684210526	4.38431782051444E-15
        0.954736842105263	4.45457561617098E-15
        1	4.479E-15
           ];

    D = interpTable(lnmo(:,1), lnmo(:,2), soc);

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
