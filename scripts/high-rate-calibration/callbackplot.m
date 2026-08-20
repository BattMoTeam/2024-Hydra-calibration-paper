function callbackplot(history, it, simulatorSetup, parameters, expdata, varargin)
%    function callbackplot(history, it, simulatorSetup, parameters, statesExp, varargin)

    opt = struct('plotEveryIt', 1, ...
                 'nonLinearSolver', [], ...
                 'objScaling', [], ...
                 'usePacked', true, ...
                 'cutoffVoltage', [], ...
                 'allowExtrapolation', true, ...
                 'doplot', true);

    opt = merge_options(opt, varargin{:});

    %vals = parameters{1}.getParameter(simulatorSetup);

    vals00 = cell(1, numel(parameters));
    for iparam = 1:numel(parameters)
        vals00{iparam} = parameters{iparam}.getParameter(simulatorSetup);
    end
    vals0 = vertcat(vals00{:});

    X = history.u{end};
    stmp = updateSetupFromScaledParameters(simulatorSetup, parameters, X);
    for iparam = 1:numel(parameters)
        vals{iparam} = parameters{iparam}.getParameter(stmp);
    end
    vals = vertcat(vals{:});

    fprintf('callbackplot it=%g\n', it);
    fprintf('vad %g\n', history.val(end));
    fprintf('u ');
    fprintf('%g ', history.u{end});
    fprintf('\n');
    fprintf('initial values ');
    fprintf('%g ', vals0);
    fprintf('\n');
    fprintf('vals ');
    fprintf('%g ', vals);
    fprintf('\n');
    fprintf('pg %g\n', history.pg(end));

    if rem(it, opt.plotEveryIt) == 0

        % Get states
        X = history.u{end};
        setup = updateSetupFromScaledParameters(simulatorSetup, parameters, X);

        if opt.usePacked
            directory = fullfile(getHydra0Dir(), 'output');
            dataFolder = md5sum(setup.model);
            problem = packSimulationProblem(setup.initstate, setup.model, setup.schedule, dataFolder, ...
                                            'Directory', directory, ...
                                            'Name', dataFolder, ...
                                            'NonLinearSolver', opt.nonLinearSolver);
            clearPackedSimulatorOutput(problem, 'Prompt', false);
            simulatePackedProblem(problem);
            [~, states] = getPackedSimulatorOutput(problem);
        else
            states = setup.run();
        end
        states = truncateStatesAtCutoff(states, opt.cutoffVoltage);

        % Quantify difference
        getTime = @(states) cellfun(@(state) state.time, states);
        getE = @(states) cellfun(@(state) state.Control.E, states);

        simulationTime = getTime(states);
        simulationVoltage = getE(states);

        if isempty(states)
            disp('RMSE unavailable: no above-cutoff simulation states');
            return
        end

        % texp = getTime(statesExp);
        % t    = getTime(states);
        % assert(norm(texp-t, 'inf') < 1e-11);

        % Eexp = getE(statesExp);
        % E    = getE(states);

        % Ediff1 = trapz(texp, abs(Eexp - E));
        % Ediff2 = sqrt(trapz(texp, (Eexp - E).^2));

        % if ~isempty(opt.objScaling)
        %     Ediff1 = Ediff1 / opt.objScaling;
        %     Ediff2 = Ediff2 / opt.objScaling;
        % end

        % str = sprintf('Integral error %g (%g)', Ediff1, Ediff2);
        % disp(str);

        if opt.allowExtrapolation
            rmse = l2error(expdata.time, expdata.U, ...
                           simulationTime, simulationVoltage, ...
                           'extrap', true);
        else
            % Use the retained simulation states as the comparison grid.
            % A very short simulation can contain several states before the
            % second experimental sample, in which case truncating the
            % experimental grid would leave fewer than two points.
            inExperimentalDomain = simulationTime >= expdata.time(1) & ...
                                   simulationTime <= expdata.time(end);
            comparisonTime = simulationTime(inExperimentalDomain);
            comparisonVoltage = simulationVoltage(inExperimentalDomain);

            if isempty(comparisonTime)
                disp('RMSE unavailable: simulation and experiment do not overlap');
                return
            end

            if isscalar(comparisonTime)
                experimentalVoltage = interp1(expdata.time, expdata.U, ...
                                              comparisonTime, 'linear');
                rmse = abs(comparisonVoltage - experimentalVoltage);
            else
                rmse = l2error(comparisonTime, comparisonVoltage, ...
                               expdata.time, expdata.U);
            end
        end
        str = sprintf('RMSE %g mV', rmse/milli);
        disp(str);

        if opt.doplot
            % Plot
            figure; hold on, grid on
            %plot(getTime(statesExp)/hour, getE(statesExp), 'displayname', 'exp')
            if opt.allowExtrapolation
                plot(expdata.time/hour, expdata.U, 'displayname', 'exp');
            else
                experimentalVoltage = interp1(expdata.time, expdata.U, ...
                                              comparisonTime, 'linear');
                plot(comparisonTime/hour, experimentalVoltage, ...
                     'displayname', 'exp');
            end
            plot(simulationTime/hour, simulationVoltage, ...
                 'displayname', 'calibrated')
            xlabel('Time  /  hour')
            ylabel('Voltage  /  V')
            legend('location', 'sw')
            title(sprintf('it=%g %s', it, str));
            drawnow
        end

    end

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
