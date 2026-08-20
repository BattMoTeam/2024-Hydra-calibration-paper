%% Script to calibrate parameters using high-rate data

clear all
close all

am    = 'ActiveMaterial';
itf   = 'Interface';
pe    = 'PositiveElectrode';
ne    = 'NegativeElectrode';
co    = 'Coating';
sd    = 'SolidDiffusion';
ctrl  = 'Control';
elyte = 'Electrolyte';
sep   = 'Separator';

mrstDebug(0);

Dnes = [1e-14, 1e-13]; % Do 1e-13 last as this is the one to use
Dnes = 1e-13;

% tag = 'no-elyte-params';
% tag = 'one-elyte-param';
% tag = 'two-elyte-params';
% tag = 'three-elyte-params';
% tag = 'elyte-bgfactor';

tags = {%'no-elyte-params', ...
        % 'one-elyte-param', ...
        % 'one-elyte-param-finers', ...
        % 'two-elyte-params', ...
        % 'three-elyte-params', ...
        'elyte-bgfactor'
       };

doplot = true;
debug = false;

for itag = 1:numel(tags)
    tag = tags{itag};
    for iDne = 1:numel(Dnes)
        Dne = Dnes(iDne);
        fprintf('Running tag=%s with Dne=%g\n', tag, Dne);

        Dpe = 1e-14;

        diary(sprintf('_diary-%s-%s-%s.txt', mfilename, tag, datestr(now, 'yyyymmdd-HHMMSS')));

        set(0, 'defaultlinelinewidth', 2)
        set(0, 'defaulttextfontsize', 15);
        set(0, 'defaultaxesfontsize', 15);


        getTime = @(states) cellfun(@(s) s.time, states);
        getE = @(states) cellfun(@(s) s.(ctrl).E, states);
        printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

        disp(tag);

        %% Fetch experimental data

        datafilename = fullfile(getHydra0Dir(), 'raw-data', 'TE_1473.mat');
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

        %commonShortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D', 'ne_bg', 'pe_bg'};
        commonShortnames = {'ne_vsa', 'pe_vsa', 'ne_D', 'pe_D'};

        switch tag
          case 'no-elyte-params'
            useRegionBruggemanCoefficients = false;
            elyteShortnames = {};

          case {'one-elyte-param', 'one-elyte-param-finer'}
            useRegionBruggemanCoefficients = false;
            elyteShortnames = {'elyte_bg'};

          case 'two-elyte-params'
            useRegionBruggemanCoefficients = true;
            elyteShortnames = {'elyte_bg_ne', 'elyte_bg_pe'};

          case 'three-elyte-params'
            useRegionBruggemanCoefficients = true;
            elyteShortnames = {'elyte_bg_ne', 'elyte_bg_pe', 'elyte_bg_sep'};

          case 'elyte-bgfactor'
            useRegionBruggemanCoefficients = true;
            elyteShortnames = {'elyte_bgfactor'};

          otherwise
            error('Unexpected tag %s', tag);
        end

        shortnames = [commonShortnames, elyteShortnames];
        disp('shortnames:');
        printer(shortnames);


        % % Estimate capacity
        % inputCap  = struct('lowRateParams'             , jsonstructEC, ...
        %                    'include_current_collectors', true);
        % outputCap = runHydra(inputCap, 'runSimulation', false);
        % cap       = computeCellCapacity(outputCap.model);

        % % Calculate Bruggeman coefficients from tortuosity and vf
        % tortuosityRef = struct(pe, 3.46, ...
        %                        ne, 3, ...
        %                        sep, 4.2);
        % bruggeman = calculateBruggemanFromTortuosity(outputCap.model, jsonstructEC, tortuosityRef);

        % if contains(tag, 'finer')
        % numTimesteps = 400;
        % else
        numTimesteps = 100;
        % end

        % printer(bruggeman)
        % keyboard;

        % Initial guess
        % input0 = struct('DRate'                         , expdata.I / cap * hour        , ...
        %                 'totalTime'                     , expdata.time(end)             , ...
        %                 'numTimesteps'                  , numTimesteps                  , ...
        %                 'lowRateParams'                 , jsonstructEC                  , ...
        %                 'Dne'                           , Dne                           , ...
        %                 'Dpe'                           , Dpe                           , ...
        %                 'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients,  ...
        %                 'include_current_collectors'    , true                          , ...
        %                 'ne_bman'                       , bruggeman.(ne).(co).bruggemanCoefficient, ...
        %                 'pe_bman'                       , bruggeman.(pe).(co).bruggemanCoefficient, ...
        %                 'sep_bman'                      , bruggeman.(sep).bruggemanCoefficient);
        input0 = struct('I'                             , expdata.I, ...
                        'totalTime'                     , expdata.time(end)             , ...
                        'numTimesteps'                  , numTimesteps                  , ...
                        'lowRateParams'                 , jsonstructEC                  , ...
                        'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients,  ...
                        'include_current_collectors'    , true);
        output0 = runHydra(input0, 'clearSimulation', false);

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
        %assert(abs(expdata.time(end) - simtimes(end)) < 1e-11);
        assert(abs(expdata.time(end) - simtimes(end))/expdata.time(end) < 3e-4);

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

        simulatorSetup = SimulationSetup(struct('model'          , output0.model   , ...
                                                'schedule'       , output0.schedule, ...
                                                'initstate'      , output0.initstate, ...
                                                'NonLinearSolver', output0.nls     , ...
                                                'OutputMinisteps', false));

        % Setup parameters to be calibrated
        HRC = HighRateCalibration(simulatorSetup, 'shortnames', shortnames);
        parameters = HRC.getParams();

        % Objective function
        lsq = @(simsetup, states, varargin) leastSquaresEI(simsetup, states, statesExp, varargin{:});
        v = lsq(simulatorSetup, output0.states);
        scaling = sum([v{:}]);
        objective = @(p, varargin) evalObjectiveBattmo(p, lsq, simulatorSetup, parameters, ...
                                                       'objScaling', scaling, varargin{:});

        if debug
            % The least squares function evaluated at the experimental values
            % should be zero
            v = lsq(simulatorSetup, statesExp);
            assert(norm([v{:}]) == 0.0);

            % Compare gradients calculated using adjoints and finite
            % difference approximation
            Xtmp = getScaledParameterVector(simulatorSetup, parameters);

            [vad, gad] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                             'gradientMethod', 'AdjointAD', ...
                                             'objScaling', scaling);

            perturbationSize = 1e-3;
            [vnum, gnum] = evalObjectiveBattmo(Xtmp, lsq, simulatorSetup, parameters, ...
                                               'gradientMethod', 'PerturbationADNUM', ...
                                               'PerturbationSize', perturbationSize, ...
                                               'objScaling', scaling);

            absErr = abs(gad - gnum);
            relErr = absErr ./ max(abs(gad) + abs(gnum), eps);
            disp(table(HRC.shortnames, gad, gnum, absErr, relErr, ...
                       'VariableNames', {'Shortname', 'Adjoint', ...
                                         'FiniteDifference', 'AbsoluteError', ...
                                         'RelativeError'}));
            fprintf('Finite-difference perturbation size in scaled coordinates: %.3e\n', ...
                    perturbationSize);
            assert(abs(vad - vnum) < eps);

            % Relative error is not meaningful for insensitive parameters,
            % so combine a relative and an absolute tolerance.
            absTol = 1e-1;
            relTol = 1e-1;
            assert(all(absErr <= absTol + relTol .* max(abs(gad), abs(gnum))), ...
                   'Adjoint and finite-difference gradients do not agree.');

            %return
        end

        %% Run optimization

        X0 = getScaledParameterVector(simulatorSetup, parameters);
        v0 = objective(X0);

        callbackfunc = @(history, it) callbackplot(history, it, simulatorSetup, parameters, expdata, ...
                                                   'plotEveryIt', 10     , ...
                                                   'objScaling' , scaling, ...
                                                   'doplot'     , doplot);

        objChangeTol = 1e-8;
        [vopt, Xopt, history] = unitBoxBFGS(X0, objective, ...
                                            'gradTol'         , 1e-3        , ...
                                            'objChangeTol'    , objChangeTol, ...
                                            'lineSearchMaxIt' , 10          , ...
                                            'maxInitialUpdate', 0.02        , ...
                                            'maximize'        , false       , ...
                                            'maxit'           , 50         , ...
                                            'logPlot'         , true        , ...
                                            'callbackfunc'    , callbackfunc, ...
                                            'plotEvolution'   , doplot);

        setupOpt = updateSetupFromScaledParameters(simulatorSetup, parameters, Xopt);

        fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, numel(history.val));
        reasonStr = getReasonStr(history);
        disp(reasonStr);

        if debug && numel(history.val) >= 2 && ...
                abs(history.val(end) - history.val(end-1)) < objChangeTol

            %% Calculate fd and adjoint gradients at final point


            [vad, gad] = evalObjectiveBattmo( ...
                Xopt, lsq, simulatorSetup, parameters, ...
                'gradientMethod', 'AdjointAD', ...
                'objScaling', scaling);

            gad = gad(:);
            shortnames = HRC.shortnames(:);

            pertsizes = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];

            nPert = numel(pertsizes);
            nParam = numel(parameters);

            vnums = nan(nPert, 1);
            gnums = cell(nPert, 1);

            for k = 1:nPert
                h = pertsizes(k);

                % Perturb in a feasible direction for each scaled parameter.
                parameterSteps = repmat({h}, nParam, 1);

                for ip = 1:nParam
                    if Xopt(ip) + h > 1
                        parameterSteps{ip} = -h;
                    end
                end

                [vnums(k), gnums{k}] = evalObjectiveBattmo( ...
                    Xopt, lsq, simulatorSetup, parameters, ...
                    'gradientMethod', 'PerturbationADNUM', ...
                    'PerturbationSize', parameterSteps, ...
                    'objScaling', scaling);

                gnums{k} = gnums{k}(:);
            end

            tbl = table(shortnames, gad, ...
                        'VariableNames', {'Shortname', 'Adjoint'});

            for k = 1:nPert
                h = pertsizes(k);
                gnum = gnums{k};

                absErr = abs(gad - gnum);

                denominator = max(abs(gad), abs(gnum));
                relErr = absErr ./ denominator;

                % Relative error is meaningless when both gradients are tiny.
                relErr(denominator < 1e-10) = NaN;

                suffix = matlab.lang.makeValidName(sprintf('h_%g', h));

                tbl.(['FD_' suffix]) = gnum;
                tbl.(['AbsErr_' suffix]) = absErr;
                tbl.(['RelErr_' suffix]) = relErr;
            end

            disp(tbl);
            keyboard;

            fprintf('\n============================================================\n');
            fprintf('BFGS STAGNATION DIAGNOSTICS\n');
            fprintf('============================================================\n');

            xPrev = history.u{end-1}(:);
            x     = Xopt(:);
            dx    = x - xPrev;

            fPrev = history.val(end-1);
            f     = history.val(end);
            df    = f - fPrev;

            %% 1. Did BFGS actually move?
            fprintf('\n[1] ACCEPTED BFGS STEP\n');

            fprintf('previous objective = %.17e\n', fPrev);
            fprintf('current objective  = %.17e\n', f);
            fprintf('objective change   = %+.17e\n', df);
            fprintf('exactly equal?     = %d\n', f == fPrev);

            fprintf('||dx||_inf         = %.17e\n', norm(dx, inf));
            fprintf('||dx||_2           = %.17e\n', norm(dx, 2));
            fprintf('identical x?       = %d\n', isequal(x, xPrev));
            fprintf('# changed params   = %d / %d\n', nnz(dx ~= 0), numel(dx));

            fprintf('line-search alpha  = %.17e\n', history.alpha(end));
            fprintf('line-search its    = %d\n', history.lsit(end));
            fprintf('line-search flag   = %d\n', history.lsfl(end));
            fprintf('history.pg         = %.17e\n', history.pg(end));

            %% 2. Recompute objective and gradient
            fprintf('\n[2] OBJECTIVE / GRADIENT AT FINAL POINT\n');

            [fRecomputed, g] = objective(x);
            g = g(:);

            fprintf('stored objective      = %.17e\n', f);
            fprintf('recomputed objective  = %.17e\n', fRecomputed);
            fprintf('difference            = %+.17e\n', fRecomputed - f);

            fprintf('||g||_inf             = %.17e\n', norm(g, inf));
            fprintf('||g||_2               = %.17e\n', norm(g, 2));


            %% 3. Which variables are on bounds?
            fprintf('\n[3] BOX CONSTRAINTS / PROJECTED GRADIENT\n');

            boundTol = 1e-10;

            atLower = x <= boundTol;
            atUpper = x >= 1-boundTol;

            % Manual projected gradient for a minimization problem.
            pg = g;

            % At lower bound, +gradient means -gradient points outside x >= 0.
            pg(atLower & g > 0) = 0;

            % At upper bound, -gradient means -gradient points outside x <= 1.
            pg(atUpper & g < 0) = 0;

            fprintf('# lower bound = %d\n', nnz(atLower));
            fprintf('# upper bound = %d\n', nnz(atUpper));
            fprintf('manual ||pg||_inf = %.17e\n', norm(pg, inf));
            fprintf('history pg         = %.17e\n', history.pg(end));

            % Physical parameter values
            setupDebug = updateSetupFromScaledParameters( ...
                simulatorSetup, parameters, x);

            physicalValues = cellfun( ...
                @(p) p.getParameter(setupDebug), parameters, ...
                'UniformOutput', false);

            physicalValues = vertcat(physicalValues{:});

            disp(table(HRC.shortnames(:), ...
                       xPrev, ...
                       x, ...
                       dx, ...
                       physicalValues, ...
                       g, ...
                       pg, ...
                       atLower, ...
                       atUpper, ...
                       'VariableNames', {'Parameter', ...
                                         'PreviousScaled', ...
                                         'CurrentScaled', ...
                                         'ScaledStep', ...
                                         'PhysicalValue', ...
                                         'Gradient', ...
                                         'ProjectedGradient', ...
                                         'AtLowerBound', ...
                                         'AtUpperBound'}));


            %% 4. Is objective evaluation repeatable?
            fprintf('\n[4] REPEATABILITY\n');

            nRepeat = 5;
            fRepeat = zeros(nRepeat, 1);

            for ir = 1:nRepeat
                fRepeat(ir) = objective(x);
            end

            disp(fRepeat);

            fprintf('repeat min      = %.17e\n', min(fRepeat));
            fprintf('repeat max      = %.17e\n', max(fRepeat));
            fprintf('repeat spread   = %.17e\n', ...
                    max(fRepeat) - min(fRepeat));

            %% 5. Objective scaling
            fprintf('\n[5] OBJECTIVE SCALING\n');

            fprintf('objScaling              = %.17e\n', scaling);
            fprintf('initial scaled obj v0   = %.17e\n', v0);
            fprintf('final scaled objective  = %.17e\n', f);
            fprintf('scaled objective change = %.17e\n', df);

            % evalObjectiveBattmo is being called with objScaling=scaling.
            % This is useful for interpreting the magnitude in original LSQ units.
            fprintf('approx raw objective change = %.17e\n', ...
                    df*scaling);


            %% 6. Objective resolution versus parameter perturbation
            fprintf('\n[6] OBJECTIVE SENSITIVITY VS STEP SIZE\n');

            hvals = [1e-1, 3e-2, 1e-2, 3e-3, ...
                     1e-3, 3e-4, 1e-4, 3e-5, ...
                     1e-5, 3e-6, 1e-6, 1e-7];

            f0debug = objective(x);

            for ip = 1:numel(x)

                fprintf('\nParameter: %s\n', HRC.shortnames{ip});
                fprintf('    h              delta-f               FD derivative\n');

                for ih = 1:numel(hvals)

                    h = hvals(ih);

                    % Pick a feasible perturbation direction.
                    if x(ip) + h <= 1
                        signedH = h;
                    elseif x(ip) - h >= 0
                        signedH = -h;
                    else
                        continue
                    end

                    xp = x;
                    xp(ip) = xp(ip) + signedH;

                    fp = objective(xp);

                    deltaF = fp - f0debug;
                    fd     = deltaF/signedH;

                    fprintf('%12.3e   %+20.12e   %+20.12e\n', ...
                            signedH, deltaF, fd);
                end
            end

            %% 7. Directional derivative test
            fprintf('\n[7] DIRECTIONAL DERIVATIVE CHECK\n');

            % Project steepest-descent direction onto feasible box directions.
            d = -g;

            d(x <= boundTol & d < 0) = 0;
            d(x >= 1-boundTol & d > 0) = 0;

            if norm(d, inf) > 0

                d = d/norm(d, inf);

                adjDirectional = g.'*d;

                fprintf('Adjoint directional derivative g''d = %.17e\n', ...
                        adjDirectional);

                hvalsDir = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];

                fprintf('    h             FD directional deriv.       ratio FD/AD\n');

                for h = hvalsDir

                    xp = x + h*d;

                    % Should already be feasible, apart from roundoff.
                    xp = min(1, max(0, xp));

                    actualH = norm(xp-x, inf);

                    if actualH == 0
                        continue
                    end

                    fp = objective(xp);

                    fdDirectional = (fp-f0debug)/h;

                    fprintf('%12.3e   %+20.12e   %+16.8e\n', ...
                            h, fdDirectional, ...
                            fdDirectional/adjDirectional);
                end

            else
                fprintf('No feasible steepest-descent direction.\n');
            end

            %% 8. Scaled -> physical parameter mapping
            fprintf('\n[8] PARAMETER MAPPING CHECK\n');

            hParam = 1e-3;

            setupBase = updateSetupFromScaledParameters( ...
                simulatorSetup, parameters, x);

            physBase = cellfun(@(p) p.getParameter(setupBase), ...
                               parameters, ...
                               'UniformOutput', false);
            physBase = vertcat(physBase{:});

            physPerturbed = nan(size(physBase));

            for ip = 1:numel(x)

                xp = x;

                if x(ip) + hParam <= 1
                    xp(ip) = x(ip) + hParam;
                else
                    xp(ip) = x(ip) - hParam;
                end

                setupP = updateSetupFromScaledParameters( ...
                    simulatorSetup, parameters, xp);

                tmp = cellfun(@(p) p.getParameter(setupP), ...
                              parameters, ...
                              'UniformOutput', false);
                tmp = vertcat(tmp{:});

                % Physical value corresponding to the parameter we moved.
                physPerturbed(ip) = tmp(ip);
            end

            physicalChange = physPerturbed - physBase;

            disp(table(HRC.shortnames(:), ...
                       x, ...
                       physBase, ...
                       physPerturbed, ...
                       physicalChange, ...
                       'VariableNames', {'Parameter', ...
                                         'ScaledValue', ...
                                         'PhysicalBase', ...
                                         'PhysicalPerturbed', ...
                                         'PhysicalChange'}));

            %% 9. Plot objective along feasible descent direction
            fprintf('\n[9] LINE SCAN\n');

            d = -pg;

            if norm(d, inf) > 0

                d = d/norm(d, inf);

                alphas = [0, logspace(-8, -1, 30)];

                fLine = nan(size(alphas));
                stepLine = nan(size(alphas));

                for ia = 1:numel(alphas)

                    xt = min(1, max(0, x + alphas(ia)*d));

                    stepLine(ia) = norm(xt-x, inf);
                    fLine(ia) = objective(xt);
                end

                figure;
                semilogx(alphas(2:end), ...
                         fLine(2:end)-fLine(1), 'o-');
                grid on
                xlabel('\alpha');
                ylabel('J(x+\alpha d)-J(x)');
                title('Objective resolution along descent direction');

                fprintf('       alpha          actual step             delta-f\n');

                for ia = 1:numel(alphas)
                    fprintf('%12.3e   %20.12e   %+20.12e\n', ...
                            alphas(ia), stepLine(ia), ...
                            fLine(ia)-fLine(1));
                end
            end
            %% 99 Adjoint and FD derivatives
            [vAdjoint, gAdjoint] = evalObjectiveBattmo(Xopt, lsq, simulatorSetup, parameters, ...
                                                       'gradientMethod', 'AdjointAD', ...
                                                       'objScaling', scaling);

            [vFDcoarse, gFDcoarse] = evalObjectiveBattmo(Xopt, lsq, simulatorSetup, parameters, ...
                                                         'gradientMethod', 'PerturbationADNUM', ...
                                                         'PerturbationSize', 1e-3, ...
                                                         'objScaling', scaling);
            [vFDfine, gFDfine] = evalObjectiveBattmo(Xopt, lsq, simulatorSetup, parameters, ...
                                                     'gradientMethod', 'PerturbationADNUM', ...
                                                     'PerturbationSize', 1e-6, ...
                                                     'objScaling', scaling);
            disp(table(HRC.shortnames, gAdjoint(:), gFDcoarse(:), gFDfine(:), ...
                       'VariableNames', {'Parameter', 'Adjoint', 'FD_1e_3', 'FD_1e_6'}));

            keyboard;
        end

        %% Extract parameters

        jsonstructHRC = HRC.export(setupOpt);
        filename = fullfile(getHydra0Dir(), 'parameters', 'high-rate-calibration-parameters.json');
        writeStruct(jsonstructHRC, filename);
        printer(jsonstructHRC);

        Dne = output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
        Dpe = output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
        filename = fullfile(getHydra0Dir(), 'parameters', sprintf('high-rate-calibration-parameters-%g-%g.json', Dne, Dpe));
        writeStruct(jsonstructHRC, filename);

        %% Run model with calibrated parameters

        % inputOpt = struct('DRate'                         , expdata.I / cap * hour        , ...
        %                   'totalTime'                     , expdata.time(end)             , ...
        %                   'numTimesteps'                  , numTimesteps                  , ...
        %                   'lowRateParams'                 , jsonstructEC                  , ...
        %                   'highRateParams'                , jsonstructHRC                 , ...
        %                   'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
        %                   'include_current_collectors'    , true);
        inputOpt = struct('I'                             , expdata.I                     , ...
                          'totalTime'                     , expdata.time(end)             , ...
                          'numTimesteps'                  , numTimesteps                  , ...
                          'lowRateParams'                 , jsonstructEC                  , ...
                          'highRateParams'                , jsonstructHRC                 , ...
                          'useRegionBruggemanCoefficients', useRegionBruggemanCoefficients, ...
                          'include_current_collectors'    , true);
        outputOpt = runHydra(inputOpt, 'clearSimulation', false);

        assert(outputOpt.model.Electrolyte.bgfactor == setupOpt.model.Electrolyte.bgfactor);

        %% Quantify differences
        vfinal = lsq(simulatorSetup, outputOpt.states);

        expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');
        tt = getTime(outputOpt.states);
        RMSE = l2error(tt, getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);

        fprintf('Final least squares values:\n');
        fprintf('vopt: %g\n', vopt);
        fprintf('Sum of squares: %g\n', sum([vfinal{:}]));
        fprintf('RMSE: %g mV\n', RMSE/milli);

        if doplot
            % plot differences
            figure; hold on; grid on;
            plot(tt, (getE(outputOpt.states) - expdataUinterp1(tt)).^2, 'displayname', '|E_{sim} - E_{exp}|^2');
            plot(tt, [vfinal{:}], 'displayname', 'vfinal');
        end

        %% Save

        Dne = output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
        Dpe = output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;

        dosavemodel = true;
        if dosavemodel
            save(sprintf('high-rate-calibrated-outputOpt-%s-%g-%g.mat', tag, Dne, Dpe));
        end

        %% Plot
        if doplot
            colors = lines(2);
            fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
            hold on;
            plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2C');
            plot(getTime(output0.states)/hour, getE(output0.states), 'color', colors(1,:), 'displayname', 'Initial guess')
            plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'Calibrated');
            xlabel('Time  /  h')
            ylabel('E  /  V')
            legend('location', 'sw')
            axis tight
            ylim([3.45, 4.9])

            dosave = true;
            if dosave
                exportgraphics(fig, sprintf('high-rate-calibration-%s-%g-%g.png', tag, Dne, Dpe), 'resolution', 300)
            end
        end

        %% Quantify difference between experiment and calibrated
        RMSE = l2error(getTime(outputOpt.states), getE(outputOpt.states), expdata.time, expdata.U, 'extrap', true);
        fprintf('wL2 error after calibration %s Dne=%g Dpe=%g: %g mV\n', tag, Dne, Dpe, RMSE/milli);

        %% Print

        fprintf('Results HRC tag=%s Dne=%g Dpe=%g\n', tag, Dne, Dpe);
        printer(jsonstructHRC);

        % Postprocess: Report effective electrode conductivities and
        % electrolyte tortuosities
        regionTags = outputOpt.model.(elyte).regionTags;
        vfs = struct();
        vfs.(ne)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 1));
        vfs.(pe)  = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 2));
        vfs.(sep) = unique(outputOpt.model.(elyte).volumeFraction(regionTags == 3));
        assert(isscalar(vfs.(ne)));
        assert(isscalar(vfs.(pe)));
        assert(isscalar(vfs.(sep)));

        if useRegionBruggemanCoefficients
            %bgfactor = outputOpt.model.(elyte).bgfactor;
            regionBruggeman = outputOpt.model.(elyte).regionBruggemanCoefficients;
            % effectiveRegionBruggeman = struct(ne , bgfactor .* regionBruggeman.(ne), ...
            %                                   pe , bgfactor .* regionBruggeman.(pe), ...
            %                                   sep, bgfactor .* regionBruggeman.(sep));
            effectiveRegionBruggeman = struct(ne , regionBruggeman.(ne), ...
                                              pe , regionBruggeman.(pe), ...
                                              sep, regionBruggeman.(sep));
       else
            bg = outputOpt.model.(elyte).bruggemanCoefficient;
            assert(isscalar(bg));
            effectiveRegionBruggeman = struct(ne, bg, pe, bg, sep, bg);
        end

        tortuosityParams = struct();
        tortuosityParams.(elyte).regionBruggemanCoefficients = effectiveRegionBruggeman;
        tau = calculateTortuosityFromBruggeman(vfs, tortuosityParams);
        disp('Effective electrolyte region Bruggeman coefficients:');
        printer(effectiveRegionBruggeman);
        disp('Tortuosities:');
        printer(tau);

        effCond = struct(pe, outputOpt.model.(pe).(co).effectiveElectronicConductivity, ...
                         ne, outputOpt.model.(ne).(co).effectiveElectronicConductivity);
        disp('Effective electronic conductivities:');
        printer(effCond);

        % For testing: print NE volumetric surface area
        fprintf('Initial diffusion Dne=%g volumetricsurfacearea=%g\n', ...
                Dne, jsonstructHRC.(ne).(co).(am).(itf).volumetricSurfaceArea);
        fprintf('RMSE %g mV\n', RMSE/milli);
        disp(reasonStr)

        diary off;

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
