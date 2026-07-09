classdef ParamOptFuncReg

    properties
        refVals % Reference solution
        I       % Current applied for the reference solution
        rawI
        relTol = 1e-10
        cap
        E
        regularization % Regularization parameters
        paramNames     % Parameter names for regularization
        currentParams  % Current parameter values for regularization
        simSetup       % Store simSetup for parameter scaling
        ScaledParams0  % Initial Params
    end

    methods

        function pof = ParamOptFuncReg(expdata, paramNames, simSetup, X0)
            pof.I = expdata.I;
            pof.rawI = expdata.rawI;
            pof.cap = expdata.cap;
            pof.E = expdata.E;
            pof.refVals = [expdata.time, expdata.E];
            pof.paramNames = paramNames;
            pof.simSetup = simSetup;
            pof = pof.initializeRegularization();
            pof.ScaledParams0 = X0;
        end

        function pof = initializeRegularization(pof)
        % Set up regularization parameters
            pof.regularization = struct();
            pof.regularization.lambda = 0.0;    % Tikhonov L2 weight

            % Remove expected values and ranges - we use scaled parameters
        end

        function obj = leastSquaresEWithReg(pof, model, states, schedule, parameters,p, varargin)

            disp('Ignoring regularization');

            % % Set current parameters first
        %     pof.setCurrentParameters(parameters);

        %     % Then call leastSquaresE with regularization enabled
        %     obj = pof.leastSquaresE(model, states, schedule, parameters,p, 'includeRegularization', true, varargin{:});

            obj = pof.leastSquaresE(model, states, schedule, varargin{:});

        end

        function setCurrentParameters(pof, params)
        % Store current parameter values for regularization
            pof.currentParams = params;
        end

        function E = evalRef(pof, t)
            E = interp1(pof.refVals(:, 1), pof.refVals(:, 2), t, 'linear', 'extrap');
            assert(isfinite(E), 'Consider allow for extrapolation in interp1');
        end

        % function obj = leastSquaresE(pof, model, states, schedule,params,pvec ,varargin)
        %     if nargin >= 5 && ~ischar(varargin{1})
        %         % First varargin is parameters
        %         parameters = varargin{1};
        %         pvec = p;

        %         varargin = varargin(2:end);
        %         pof.currentParams = pvec;  % Store for regularization
        %     end

        %     opt = struct('ComputePartials', false, ...
        %                  'tStep'          , []   , ...
        %                  'state'          , []   , ...
        %                  'from_states'    , true, ...
        %                  'includeRegularization', false);
        %     opt = merge_options(opt, varargin{:});

        %     dts = schedule.step.val;

        %     if isempty(opt.tStep) %do all
        %         time     = cumsum(dts);
        %         numSteps = numel(dts);
        %     else
        %         time     = sum(dts(1:(opt.tStep)));
        %         numSteps = 1;
        %         dts      = dts(opt.tStep);
        %     end

        %     obj = repmat({[]}, numSteps, 1);

        %     for i = 1:length(states)
        %         II(i) = states{i}.Control.I;
        %     end
        %     capn = abs(trapz(time, II));

        %     % Calculate total regularization once using SCALED parameters
        %     total_reg_obj = 0;
        %     initScaledParams = pof.ScaledParams0;
        %     if opt.includeRegularization
        %         % Get scaled parameters for consistent regularization
        %         scaledParams =  pvec;
        %         total_reg_obj = pof.computeRegularization(scaledParams, initScaledParams);
        %     end

        %     % Distribute regularization evenly across time steps
        %     reg_obj_per_step = total_reg_obj / numSteps;

        %     w_Ee = 0;
        %     w_Ie = 0;
        %     w_ce = 0;
        %     totTime =0;

        %     for k = 1 : min(length(states),numSteps)
        %         t  = time(k);
        %         dt = dts(k);
        %         w_Ee = w_Ee + pof.evalRef(t);
        %         w_Ie = w_Ie + pof.rawI(k);
        %         w_ce = w_ce + pof.cap;
        %         totTime = totTime +dt;
        %     end
        %     w_E = 1/w_Ee;
        %     w_I = 1/w_Ie;
        %     w_c = 1/w_ce;
        %     for k = 1 : numSteps
        %         t  = time(k);
        %         dt = dts(k);

        %         % Find states for time t (given by schedule)
        %         state = pof.findState(t, states, dt);

        %         if opt.ComputePartials
        %             if opt.from_states
        %                 state = model.initStateAD(state);
        %             else
        %                 state = opt.state;
        %             end
        %         end

        %         Ee = state.Control.E;
        %         Ic = state.Control.I;

        %         Eref = pof.evalRef(t);
        %         Iref = pof.rawI(k);
        %         capref = pof.cap;

        %         % Data misfit for this time step
        %         %data_obj = dt*(w_E*((Ee - Eref))^2 + 0.1.*w_I*((Ic-Iref))^2 + w_c*((capref-capn))^2);
        %         data_obj = (Ee - Eref)^2./Ee^2 * dt + 0.001.*(Ic-Iref)^2./Ic^2* dt +(capref-capn)^2./capn^2./dt;
        %         % Add regularization portion to EVERY time step
        %         obj{k} = data_obj + reg_obj_per_step;
        %     end
        % end

        function obj = leastSquaresE(pof, model, states, schedule, varargin)

            opt = struct('ComputePartials', false, ...
                         'tStep'          , []   , ...
                         'state'          , []   , ...
                         'from_states'    , true);
            opt = merge_options(opt, varargin{:});

            dts = schedule.step.val;

            if isempty(opt.tStep) %do all
                time     = cumsum(dts);
                numSteps = numel(dts);
            else
                time     = sum(dts(1:(opt.tStep)));
                numSteps = 1;
                dts      = dts(opt.tStep);
            end

            obj = repmat({[]}, numSteps, 1);

            for k = 1 : numSteps

                t  = time(k);
                dt = dts(k);

                % Find states for time t (given by schedule)
                state = pof.findState(t, states, dt);

                if opt.ComputePartials
                    if opt.from_states
                        state = model.initStateAD(state);
                    else
                        state = opt.state;
                    end
                end

                E = state.Control.E;

                Eref = pof.evalRef(t);

                obj{k} = (E - Eref)^2 * dt;

            end

        end

        function reg_obj = computeRegularization(pof, scaledParameters, scaledParameters0)
        % Compute Tikhonov L2 regularization using SCALED parameters
            reg_obj = 0;

            if isempty(scaledParameters) || isempty(pof.regularization)
                return;
            end

            % Use scaled parameters directly (they're already normalized)
            scaledValues = scaledParameters;

            scaledValues0 = scaledParameters0;
            if isa(scaledValues, 'ADI')
                scaledValues = value(scaledValues);
            end

            % Simple Tikhonov L2: λ * ||u||² where u are scaled parameters
            lambda = pof.regularization.lambda;
            if iscell(scaledValues)
                for i = 1:length(scaledValues)
                    scaledVal = scaledValues{1}(i)-scaledValues0(i);
                    reg_obj = reg_obj + lambda * scaledVal^2;
                end
            else
                for i = 1:length(scaledValues)
                    scaledVal = scaledValues(i)-scaledValues0(i);
                    reg_obj = reg_obj + lambda * scaledVal^2;

                end
            end


            % Optional: Add penalty for extreme values (outside [-2, 2] in scaled space)
            % for i = 1:length(scaledValues)
            %     scaledVal = scaledValues(i);
            %     if abs(scaledVal) > 2.0
            %         reg_obj = reg_obj + lambda * 10.0 * (abs(scaledVal) - 2.0)^2;
            %     end
            % end
        end

        function currentValues = extractParamsSafely(pof, model)
        % Safe parameter extraction with warning suppression
            warning('off', 'ADI:doubleDeprecated');

            try
                % Pass the model to getfun
                currentValues = pof.currentParams{1}.getfun(model, []);
                if isa(currentValues, 'ADI')
                    currentValues = value(currentValues);
                end
            catch
                % Fallback: return zeros
                currentValues = zeros(1, length(pof.paramNames));
            end

            warning('on', 'ADI:doubleDeprecated');
        end

        function state = findState(pof, tstar, states, dt)
            states0 = states;
            idx = cellfun(@(x) isempty(x), states);
            states = states(~idx);
            times = cellfun(@(x) x.time, states);
            [tdiff, idx] = min(abs(times - tstar));
            assert(tdiff / dt < pof.relTol);
            state = states{idx};
        end

        function setRegularizationWeight(pof, lambda)
            pof.regularization.lambda = lambda;
        end

        function setSimSetup(pof, simSetup)
            pof.simSetup = simSetup;
        end

    end

end
