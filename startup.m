function startup()

    repoRoot = fileparts(mfilename('fullpath'));
    dirnames = {'raw-data', 'scripts', 'parameters', 'utils'};

    for idir = 1:numel(dirnames)
        addpath(genpath(fullfile(repoRoot, dirnames{idir})));
    end

    cwdir = pwd();
    battmo = resolveBattMoPath(repoRoot);

    cleanupObj = onCleanup(@() cd(cwdir)); %#ok<NASGU>
    configurePythonExecutable(repoRoot);
    cd(battmo);
    run('startupBattMo.m')
    cd(cwdir);

    mrstModule add ad-core optimization mpfa


    fprintf('\nCurrent directory: %s\n\n', pwd());

end


function battmo = resolveBattMoPath(repoRoot)

    parent1 = fileparts(repoRoot);
    parent2 = fileparts(parent1);

    candidates = {};

    % Preferred path for CI/users: set this explicitly to your BattMo clone.
    envBattmo = getenv('BATTMO_DIR');
    if ~isempty(envBattmo)
        candidates{end+1} = envBattmo; %#ok<AGROW>
    end

    % Common workspace layouts used in this repository.
    candidates{end+1} = fullfile(parent2, 'BattMo', 'BattMo', 'BattMo'); %#ok<AGROW>
    candidates{end+1} = fullfile(parent2, 'BattMo', 'BattMo'); %#ok<AGROW>
    candidates{end+1} = fullfile(parent2, 'BattMo'); %#ok<AGROW>
    candidates{end+1} = fullfile(parent1, 'BattMo'); %#ok<AGROW>
    candidates{end+1} = fullfile(repoRoot, 'BattMo'); %#ok<AGROW>

    candidates = unique(candidates, 'stable');

    for i = 1:numel(candidates)
        c = candidates{i};
        if isfolder(c) && isfile(fullfile(c, 'startupBattMo.m'))
            battmo = c;
            fprintf('Using BattMo path: %s\n', battmo);
            return;
        end
    end

    error(['Could not find BattMo startup script. Set BATTMO_DIR to the ', ...
           'BattMo root containing startupBattMo.m. Checked candidates:\n%s'], ...
          strjoin(candidates, '\n'));

end


function configurePythonExecutable(repoRoot)

    pe = pyenv();
    if pe.Version ~= ""
        fprintf('Using Python executable: %s\n', pe.Executable);
        return;
    end

    candidates = {};

    envvars = {'HYDRA_PYTHON_EXECUTABLE', ...
               'PYTHON_EXECUTABLE'      , ...
               'BATTMO_PYTHON_EXECUTABLE'};

    for i = 1:numel(envvars)
        envval = getenv(envvars{i});
        if ~isempty(envval)
            candidates{end+1} = envval; %#ok<AGROW>
        end
    end

    if ispc
        localPrograms = fullfile(getenv('LOCALAPPDATA'), 'Programs', 'Python');
        candidates = [candidates, ...
                      {fullfile(repoRoot, 'env', 'Scripts', 'python.exe')} ...
                      {fullfile(localPrograms, 'Python313', 'python.exe')} ...
                      {fullfile(localPrograms, 'Python312', 'python.exe')} ...
                      {fullfile(localPrograms, 'Python311', 'python.exe')} ...
                      {fullfile(localPrograms, 'Python310', 'python.exe')} ...
                      {'C:\Program Files\Python313\python.exe'} ...
                      {'C:\Program Files\Python312\python.exe'} ...
                      {'C:\Program Files\Python311\python.exe'} ...
                      {'C:\Program Files\Python310\python.exe'} ...
                      {'C:\Program Files\Python39\python.exe'} ...
                      {'C:\Program Files\Python38\python.exe'}];
    else
        candidates = [candidates, ...
                      {fullfile(repoRoot, 'env', 'bin', 'python')} ...
                      {'/usr/bin/python3.13'} ...
                      {'/usr/bin/python3.12'} ...
                      {'/usr/bin/python3.11'} ...
                      {'/usr/bin/python3.10'} ...
                      {'/usr/bin/python3.9'} ...
                      {'/usr/bin/python3.8'}];
    end

    candidates = unique(candidates, 'stable');

    for i = 1:numel(candidates)
        candidate = candidates{i};
        if ~isfile(candidate)
            continue;
        end

        try
            pyenv('Version', candidate);
            fprintf('Using Python executable: %s\n', pyenv().Executable);
            return;
        catch
            % Keep trying until we find a Python supported by this MATLAB release.
        end
    end

    warning(['No supported Python executable configured. JSON validation in BattMo ', ...
             'will fail until pyenv is set explicitly. Set HYDRA_PYTHON_EXECUTABLE ', ...
             'or call pyenv(''Version'', <python-executable>) before running scripts.']);

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
