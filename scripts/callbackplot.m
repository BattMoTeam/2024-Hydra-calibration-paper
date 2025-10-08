function callbackplot(history, it, simulatorSetup, parameters, statesExp, varargin)

    opt = struct('plotEveryIt', 1, ...
                 'nonLinearSolver', [], ...
                 'objScaling', []);

    opt = merge_options(opt, varargin{:});

    vals = parameters{1}.getParameter(simulatorSetup);

    fprintf('callbackplot it=%g\n', it);
    fprintf('vad %g\n', history.val(end));
    fprintf('u ');
    fprintf('%g ', history.u{end});
    fprintf('\n');
    fprintf('vals ');
    fprintf('%g ', vals);
    fprintf('\n');
    fprintf('pg %g\n', history.pg(end));

    if rem(it, opt.plotEveryIt) == 0

        % Get states
        X = history.u{end};
        setup = updateSetupFromScaledParameters(simulatorSetup, parameters, X);

        dataFolder = md5sum(setup.model);
        problem = packSimulationProblem(setup.state0, setup.model, setup.schedule, dataFolder, ...
                                           'NonLinearSolver', opt.nonLinearSolver);
        clearPackedSimulatorOutput(problem, 'Prompt', false);
        simulatePackedProblem(problem);
        [~, states] = getPackedSimulatorOutput(problem);

        % Quantify difference
        getTime = @(states) cellfun(@(state) state.time, states);
        getE = @(states) cellfun(@(state) state.Control.E, states);

        texp = getTime(statesExp);
        t    = getTime(states);
        assert(norm(texp-t, 'inf') < 1e-12);

        Eexp = getE(statesExp);
        E    = getE(states);

        Ediff1 = trapz(texp, abs(Eexp - E));
        Ediff2 = trapz(texp, (Eexp - E).^2);

        if ~isempty(opt.objScaling)
            Ediff1 = Ediff1 / opt.objScaling;
            Ediff2 = Ediff2 / opt.objScaling;
        end

        str = sprintf('Integral error %g (%g)', Ediff1, Ediff2);
        disp(str);

        % Plot
        figure; hold on, grid on
        plot(getTime(statesExp)/hour, getE(statesExp), 'displayname', 'exp')
        plot(getTime(states)/hour, getE(states), 'displayname', 'calibrated')
        xlabel('time / hour')
        ylabel('voltage / V')
        legend('location', 'sw')
        title(sprintf('it=%g %s', it, str));
        drawnow

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
