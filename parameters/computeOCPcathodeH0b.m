function [OCP, dUdT] = computeOCPcathodeH0b(c, ~, cmax)

    lnmo_ds = [
        0.140000000000000   4.875695054478530
        0.146026026026026   4.827470014308610
        0.155495495495496   4.797539228965480
        0.162382382382382   4.785999205191700
        0.170990990990991   4.776375216330410
        0.195095095095095   4.760944464890990
        0.212312312312312   4.755003870000000
        0.238998998998999   4.750243314467020
        0.363823823823824   4.742038624274030
        0.436136136136136   4.732696020736200
        0.473153153153153   4.723104856656520
        0.511031031031031   4.703714883231080
        0.525665665665666   4.699329579529650
        0.647907907907908   4.687306302739120
        0.699559559559560   4.677593726400820
        0.733993993993994   4.665391142188140
        0.765845845845846   4.646515269734150
        0.789949949949950   4.622604094161990
        0.808028028028028   4.592155089701570
        0.819219219219219   4.559066248141210
        0.826106106106106   4.520866785419220
        0.830410410410410   4.469804330354880
        0.836436436436437   4.325886328703930
        0.840740740740741   4.293299192952090
        0.847627627627628   4.261850621446690
        0.866566566566567   4.201906063397570
        0.906166166166166   4.097083219658340
        0.925965965965966   4.053347058574070
        0.962122122122122   3.982634173975450
        0.975895895895896   3.950633224366950
        0.985365365365365   3.919831925017250
        0.991391391391391   3.885542787977680
        0.993973973973974   3.853994257271790
        0.995695695695696   3.781576034887230
        1.000000000000000   2.995543758834360
              ];

    stoc = c/cmax;
    OCP = interpTable(lnmo_ds(:,1), lnmo_ds(:,2), stoc);
    dUdT = 0;

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
