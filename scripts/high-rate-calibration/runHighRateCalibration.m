% Script to calibrate parameters using high-rate data

clear all
close all

neDs = 1e-13; %[1e-13, 1e-14];
% tag = 'no-elyte-params';
% tag = 'one-elyte-param';
% tag = 'two-elyte-params';
% tag = 'three-elyte-params'; % 1e-13 goes to it=150

tags = {%'no-elyte-params', ...
        % 'one-elyte-param', ...
        % %'one-elyte-param-finers', ...
        % 'two-elyte-params', ...
        'three-elyte-params'
       };

doplot = false;

for itag = 1:numel(tags)
    tag = tags{itag};
    for ineD = 1:numel(neDs)
        neD = neDs(ineD);
        fprintf('Running tag=%s with neD=%g\n', tag, neD);

        peD = 1e-14;

        diary(sprintf('_diary-%s-%s-%s.txt', mfilename, tag, datestr(now, 'yyyymmdd-HHMMSS')));

        mrstDebug(0);

        set(0, 'defaultlinelinewidth', 2)
        set(0, 'defaulttextfontsize', 15);
        set(0, 'defaultaxesfontsize', 15);

        am    = 'ActiveMaterial';
        itf   = 'Interface';
        pe    = 'PositiveElectrode';
        ne    = 'NegativeElectrode';
        co    = 'Coating';
        sd    = 'SolidDiffusion';
        ctrl  = 'Control';
        elyte = 'Electrolyte';
        sep   = 'Separator';

        getTime = @(states) cellfun(@(s) s.time, states);
        getE = @(states) cellfun(@(s) s.(ctrl).E, states);
        printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

        debug = false;

        % tag = 'no-elyte-params';
        % tag = 'one-elyte-param';
        % tag = 'two-elyte-params';
        % tag = 'three-elyte-params'; % 1e-13 goes to it=150
        disp(tag);

        %% Fetch experimental data

        % Original data
        datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
        saveddata    = load(datafilename);
        dataraw      = saveddata.experiment;

        % Highest DRate is last
        k = numel(dataraw.time);
        expdata = struct('time', dataraw.time{k} * hour, ...
                         'U'   , dataraw.voltage{k}    , ...
                         'I'   , abs(mean(dataraw.current{k})));

        %% Initial guess using equilibrium calibration data

        filename     = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
        jsonstructEC = parseBattmoJson(filename);

        switch tag
          case {'no-elyte-params', 'one-elyte-param', 'one-elyte-param-finer'}
            useRegionBruggemanCoefficients = false;
          case {'two-elyte-params', 'three-elyte-params'}
            useRegionBruggemanCoefficients = true;
          otherwise
            error('Unexpected tag %s', tag);
        end

        % Find capacity
        inputCap  = struct('lowRateParams'             , jsonstructEC, ...
                           'include_current_collectors', true);
        outputCap = runHydra(inputCap, 'runSimulation', false);
        cap       = computeCellCapacity(outputCap.model);

        % Diffusion must be a scalar when calibrated
        % neD = mean(computeDanodeH0b(linspace(0, 1, 100)));
        % peD = mean(computeDcathodeH0b(linspace(0.14, 1, 100)));

        % Both electrodes with Ds = 1e-14
        % Graphite: Ds= 1e-13 and LNMO: Ds=1e-14
        % neD = 1e-14;
        % peD = 1e-14;
        % neD = [];
        % peD = [];

        % Calculate Bruggeman coefficients from tortuosity and vf
        tortuosityRef = struct(pe, 3.46, ...
                               ne, 3, ...
                               sep, 4.2);
        bruggeman = calculateBruggemanFromTortuosity(outputCap.model, jsonstructEC, tortuosityRef);

        if contains(tag, 'finer')
            numTimesteps = 400;
        else
            numTimesteps = 100;
        end

        % Initial guess
        input0 = struct('DRate'                         , expdata.I / cap * hour        , ...
                        'totalTime'                     , expdata.time(end)             , ...
                        'numTimesteps'                  , numTimesteps                  , ...
                        'lowRateParams'                 , jsonstructEC                  , ...
                        'neD'                           , neD                           , ...
                        'peD'                           , peD                           , ...
                        'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients,  ...
                        'include_current_collectors'    , true                          , ...
                        'ne_bman'                       , bruggeman.(ne)                , ...
                        'pe_bman'                       , bruggeman.(pe)                , ...
                        'sep_bman'                      , bruggeman.(sep));
        output0 = runHydra(input0, 'clearSimulation', true);

        if debug
            % Check how exp and initial guess compare
            figure; hold on; grid on;
            plot(expdata.time/hour, expdata.U, 'k--');
            plot(getTime(output0.states)/hour, getE(output0.states));
            xlabel('time / h')
            ylabel('potential / V')
            title('initial guess')
            drawnow
        end

        %% Setup optimization

        % Evaluate experimental data at simulation times (allow for
        % extrapolation since expdata.time(end) is very close to
        % output.states{end}.time)
        simtimes = getTime(output0.states);
        assert(expdata.time(1) <= simtimes(1));
        assert(abs(expdata.time(end) - simtimes(end)) < 1e-11);

        Evals     = interp1(expdata.time, expdata.U, simtimes, 'linear', 'extrap');
        statesExp = cell(numel(output0.states), 1);

        for k = 1:numel(output0.states)
            statesExp{k}.time     = simtimes(k);
            statesExp{k}.(ctrl).E = Evals(k);
        end

        if debug
            % Check that the extracted values are the same as the raw values
            figure; hold on; grid on;
            plot(expdata.time/hour, expdata.U, 'k--');
            plot(getTime(statesExp)/hour, getE(statesExp));
            xlabel('Time / h')
            ylabel('Potential / V')
            title('statesExp')
            drawnow
        end

        simulatorSetup = struct('model'   , output0.model   , ...
                                'schedule', output0.schedule, ...
                                'state0'  , output0.initstate);

        % Setup parameters to be calibrated
        HRC = HighRateCalibration(simulatorSetup, tag);
        parameters = HRC.getParams();

        % Objective function
        lsq = @(model, states, schedule, varargin) leastSquaresEI(model, states, statesExp, schedule, varargin{:});
        v = lsq(simulatorSetup.model, output0.states, simulatorSetup.schedule);
        scaling = sum([v{:}]);

        objective = @(p, varargin) evalObjectiveBattmo(p, lsq, simulatorSetup, parameters, ...
                                                       'objScaling', scaling, varargin{:});

        if debug
            % The least squares function evaluated at the experimental values
            % should be zero
            v = lsq(output0.model, statesExp, simulatorSetup.schedule);
            assert(norm([v{:}]) == 0.0);

            % Compare gradients calculated using adjoints and finite
            % difference approximation
            Xtmp = getScaledParameterVector(simulatorSetup, parameters);

            [vad, gad] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                             'gradientMethod', 'AdjointAD', ...
                                             'objScaling', scaling);

            [vnum, gnum] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                               'gradientMethod', 'PerturbationADNUM', ...
                                               'PerturbationSize', 1e-7, ...
                                               'objScaling', scaling);
            assert(abs(vad - vnum) < eps);
            assert(all(abs(gad) > 0));
            assert(all(abs(gnum) > 0));
            assert(norm((gad-gnum)./gnum, 'inf') < 1e-3);

        end

        %% Run optimization

        X0 = getScaledParameterVector(simulatorSetup, parameters);
        v0 = objective(X0);

        callbackfunc = @(history, it) callbackplot(history, it, simulatorSetup, parameters, statesExp, ...
                                                   'plotEveryIt', 10, ...
                                                   'objScaling', scaling, ...
                                                   'doplot', doplot);

        [vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                            'objChangeTol', 1e-8 , ...
                                            'maximize'    , false, ...
                                            'maxit'       , 150  , ...
                                            'logPlot'     , true, ...
                                            'callbackfunc', callbackfunc, ...
                                            'plotEvolution', doplot);

        % [vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
        %                                     'objChangeTol', 1e-14 , ...
        %                                     'gradTol'     , 1e-5  , ...
        %                                     'maximize'    , false, ...
        %                                     'maxit'       , 150  , ...
        %                                     'logPlot'     , true, ...
        %                                     'callbackfunc', callbackfunc);

        setupOpt = updateSetupFromScaledParameters(simulatorSetup, parameters, Xopt);

        fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));
        reasonStr = getReasonStr(history);
        disp(reasonStr);

        %% Extract parameters

        jsonstructHRC = HRC.export(setupOpt);
        filename = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
        writeJsonStruct(jsonstructHRC, filename);
        printer(jsonstructHRC);

        neD = output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
        peD = output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
        filename = fullfile(getHydra0Dir(), 'parameters', sprintf('high-rate-calibration-parameters-%g-%g.json', neD, peD));
        writeJsonStruct(jsonstructHRC, filename);

        %% Run model with calibrated parameters


        inputOpt = struct('DRate'                         , expdata.I / cap * hour        , ...
                          'totalTime'                     , expdata.time(end)             , ...
                          'numTimesteps'                  , numTimesteps                  , ...
                          'lowRateParams'                 , jsonstructEC                  , ...
                          'highRateParams'                , jsonstructHRC                 , ...
                          'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                          'include_current_collectors'    , true);
        outputOpt = runHydra(inputOpt);

        %% Quantify differences
        vfinal = lsq(outputOpt.model, outputOpt.states, simulatorSetup.schedule);

        expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');
        tt = getTime(outputOpt.states);
        L2 = trapz(tt, (getE(outputOpt.states) - expdataUinterp1(tt)).^2) / tt(end);

        fprintf('Final least squares values:\n');
        fprintf('vopt: %g\n', vopt);
        fprintf('Sum of squares: %g\n', sum([vfinal{:}]));
        fprintf('L2^2 error: %g\n', L2);

        if doplot
            % plot differences
            figure; hold on; grid on;
            plot(tt, (getE(outputOpt.states) - expdataUinterp1(tt)).^2, 'displayname', '|E_{sim} - E_{exp}|^2');
            plot(tt, [vfinal{:}], 'displayname', 'vfinal');
        end

        %% Save

        neD = output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
        peD = output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;

        dosavemodel = true;
        if dosavemodel
            save(sprintf('high-rate-calibrated-outputOpt-%s-%g-%g.mat', tag, neD, peD));
        end

        %% Plot
        if doplot
            colors = lines(2);
            fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
            hold on;
            plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2 C');
            plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'Initial guess')
            plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'Calibrated');
            xlabel('Time  /  h')
            ylabel('E  /  V')
            legend('location', 'sw')
            axis tight
            ylim([3.45, 4.9])

            dosave = true;
            if dosave
                exportgraphics(fig, sprintf('high-rate-calibration-%s-%g-%g.png', tag, neD, peD), 'resolution', 300)
            end
        end

        %% Quantify difference between experiment and calibrated
        tt = getTime(outputOpt.states);
        wL2 = sqrt(trapz(tt, (getE(outputOpt.states) - expdataUinterp1(tt)).^2)) / tt(end) / milli;
        fprintf('wL2 error after calibration %s Dne=%g Dpe=%g: %g mV\n', tag, neD, peD, wL2);


        %% Print

        fprintf('Results HRC tag=%s neD=%g peD=%g\n', tag, neD, peD);
        printer(jsonstructHRC);

        % Postprocess: Report effective electrode conductivities and
        % electrolyte tortuosities
        tau = calculateTortuosityFromBruggeman(outputOpt.model, jsonstructHRC);
        disp('Tortuosities:');
        printer(tau);

        effCond = struct(pe, outputOpt.model.(pe).(co).effectiveElectronicConductivity, ...
                         ne, outputOpt.model.(ne).(co).effectiveElectronicConductivity);
        disp('Effective electronic conductivities:');
        printer(effCond);

        diary off;

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

    end
end
