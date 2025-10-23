function startup()

    dirnames = {'rawData', 'scripts', 'parameters', 'utils'};

    for idir = 1:numel(dirnames)
        addpath(genpath(dirnames{idir}));
    end

    username = getUsername();
    cwdir = pwd();

    switch username

      case 'august'

        battmo = '/home/august/Projects/Battery/BattMo';

      case 'xavier'

        battmo = '/home/xavier/Matlab/Projects/battmo';

      otherwise

        error('Unknown user name %s', username);

    end

    cd(battmo);
    run('startupBattMo.m')
    cd(cwdir);

    mrstModule add ad-core optimization mpfa


    fprintf('\nCurrent directory: %s\n\n', pwd());

end


function status = syscall(str)

    status = system(str);
    assert(status == 0, str);

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
