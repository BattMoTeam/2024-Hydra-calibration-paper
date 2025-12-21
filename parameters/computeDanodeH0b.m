function D = computeDanodeH0b(soc)

    graphite = [0   1.92E-13
                0.1 1.16E-13
                0.2 7.12E-14
                0.3 4.50E-14
                0.4 3.00E-14
                0.5 2.20E-14
                0.6 1.82E-14
                0.7 1.73E-14
                0.8 1.85E-14
                0.9 2.13E-14
                1   2.57E-14
               ];

    D = interpTable(graphite(:,1), graphite(:,2), soc);

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
