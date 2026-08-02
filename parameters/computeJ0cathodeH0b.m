function j0 = computeJ0cathodeH0b(soc)

    lnmo = [
        0.000000000000000   0.000146947187111
        0.025000000000000   0.000143560314165
        0.049999999999999   0.000096574397799
        0.099999999999999   0.000033828337612
        0.149999999999999   0.000262322367243
        0.199999999999999   0.000395526276763
        0.249999999999999   0.000439698519876
        0.299999999999999   0.000446259692946
        0.350000000000000   0.000432657905944
        0.400000000000000   0.000415872886922
        0.450000000000000   0.000392415831529
        0.500000000000000   0.000367849409319
        0.550000000000000   0.000391915976034
        0.600000000000000   0.000395092933335
        0.650000000000000   0.000390933653701
        0.700000000000000   0.000380879000213
        0.750000000000000   0.000377009511009
        0.800000000000000   0.000359439847188
        0.850000000000000   0.000344566811475
        0.900000000000000   0.000315637394329
        0.950000000000000   0.000286408292036
        0.975000000000000   0.000210757234065
        1.000000000000000   0.000403310274857
           ];

    j0 = interpTable(lnmo(:,1), lnmo(:,2)/centi^2, soc);

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
