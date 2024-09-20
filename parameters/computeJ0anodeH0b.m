function j0 = computeJ0anodeH0b(soc)

    graphite = [0	6.14180704504965E-05
                0.025	0.000195048879542242
                0.05	0.000285609207829683
                0.1	0.000287777671775122
                0.15	0.000357833797055224
                0.2	0.000375161042355194
                0.25	0.000399042436699301
                0.3	0.000385888124406325
                0.35	0.000391077216469399
                0.4	0.000396434407803536
                0.45	0.000402950087744402
                0.5	0.000419105848461824
                0.55	0.000504918252889569
                0.6	0.000470739079700868
                0.65	0.000468047272216724
                0.7	0.000467314920246704
                0.75	0.000471058023371167
                0.8	0.000476839651936555
                0.85	0.000481773917337519
                0.9	0.000486312366426169
                0.95	0.000497030884387319
                0.975	0.000506113895934315
                1	0.000552665137750124
               ];

    j0 = interpTable(graphite(:,1), graphite(:,2)/centi^2, soc);

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
