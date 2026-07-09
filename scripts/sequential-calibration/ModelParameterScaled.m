classdef ModelParameterScaled
    properties
        name
        type          = 'value';    % 'value'/'multiplier'
        boxLims                     % upper/lower value(s) for parameters (used for scaling)
        subset                      % subset of parameters (or subset of wells)
        scaling       = 'linear'    % 'linear'/'log'/'exp'/'custom'
        referenceValue              % parameter reference values (used for type 'multiplier') 
        belongsTo                   % model/shcedule/state0/
        location                    % e.g., {'model', 'operators', 'T'}
        nParam                      % number of parameters
        lumping                     % parameter lumping vector (partition vector) 
        setfun                      % possibly custom set-function (default is setfield)
        getfun                      % (default getfield)
        scalingBase   = nan;          
        controlSteps  = [];         % Default set/get for all control steps
        controlType   = 'none';     % types 'bhp', 'rate', 'orat', ... and 'policy' requires special 
                                    % treatment for sensitivitry computations
        wellSubsets   = [];         % subset of individual well parameters 
        extraLocations = [];        % additional locations of parameter
        
        % Enhanced properties for individual parameter scaling
        scalingFunction = [];       % Custom scaling function handle
        unscalingFunction = [];     % Custom unscaling function handle
        gradientScalingFunction = []; % Custom gradient scaling function handle
        parameterType = 'generic';  % Type identifier for parameter-specific handling
        
        % NEW: Individual scaling for each parameter
        individualScaling = [];     % Cell array of scaling types for each parameter
        individualBoxLims = [];     % Individual box limits for each parameter
        parameterNames = {};        % Names of individual parameters (for 7 parameters)
        shortnames = {}
    end
    
    methods
        function p = ModelParameterScaled(setup, varargin)
            [p, extra] = merge_options(p, varargin{:});
            assert(~isempty(p.name), 'Parameter name can''t be defaulted');
            
            % Set parameter type based on name pattern
            p = setParameterType(p);
            
            if isempty(p.belongsTo) || isempty(p.location) 
                p = setupByName(p, setup);
            end
            if isempty(p.setfun)% use default
                p.setfun = @setfield;
            end
            if isempty(p.getfun)% use default
                p.getfun = @getfield;
            end
            
            % Setup custom scaling functions based on parameter type
            p = setupCustomScaling(p);
            
            opt = struct('relativeLimits', [.5 2], ...
                         'uniformLimits',  true,   ...
                         'wellSubsets',     []);
            opt = merge_options(opt, extra{:});
            p   = setupDefaults(p, setup, opt);
            
            % NEW: Setup individual parameter scaling for your 7 parameters
            p = setupIndividualParameterScaling(p);
            
            checkSetup(p, setup);
        end
        
        %------------------------------------------------------------------
        function vs = scale(p, pval)
            % map parameter pval to "control"-vector v \in [0,1]
            if ~isempty(p.individualScaling)
                % NEW: Scale each parameter individually
                vs = zeros(size(pval));
                for i = 1:min(numel(pval), p.nParam)
                    single_pval = pval(i);
                    if i <= numel(p.individualScaling)
                        scaling_type = p.individualScaling{i};
                        box_lims = p.individualBoxLims(i,:);
                    else
                        scaling_type = p.scaling;
                        box_lims = p.boxLims(min(i, size(p.boxLims,1)),:);
                    end
                    vs(i) = p.scaleSingleParameter(single_pval, scaling_type, box_lims,i);
                end
            elseif ~isempty(p.scalingFunction)
                % Use custom scaling function
                vs = p.scalingFunction(pval, p);
            else
                % Use built-in scaling
                vs = (pval-p.boxLims(:,1))./diff(p.boxLims, [], 2) +eps;
                if strcmp(p.scaling, 'exp')
                    vs = expScale(vs, p,[]);
                elseif strcmp(p.scaling, 'log')
                    vs = logScale(vs, p,[]);
                end
            end
        end
        
        %------------------------------------------------------------------
        function pval = unscale(p, vs)
            % retrieve parameter pval from "control"-vector v \in [0,1]
            if ~isempty(p.individualScaling)
                % NEW: Unscale each parameter individually
                pval = zeros(size(vs));
                for i = 1:min(numel(vs), p.nParam)
                    single_vs = vs(i);
                    if i <= numel(p.individualScaling)
                        scaling_type = p.individualScaling{i};
                        box_lims = p.individualBoxLims(i,:);
                    else
                        scaling_type = p.scaling;
                        box_lims = p.boxLims(min(i, size(p.boxLims,1)),:);
                    end
                    pval(i) = p.unscaleSingleParameter(single_vs, scaling_type, box_lims, i);
                end
            elseif ~isempty(p.unscalingFunction)
                % Use custom unscaling function
                pval = p.unscalingFunction(vs, p);
            else
                % Use built-in unscaling
                pval = vs;
                if strcmp(p.scaling, 'exp')
                    pval = logScale(pval, p,[]);
                elseif strcmp(p.scaling, 'log')
                    pval = expScale(pval, p,[]);
                end
                pval = pval.*diff(p.boxLims, [], 2) + p.boxLims(:,1);
            end
        end
        
        %------------------------------------------------------------------
        function gs = scaleGradient(p, g, pval)
            % map gradient wrt param to gradient vs "control"-vector
            if ~isempty(p.individualScaling)
                % NEW: Scale gradient for each parameter individually
                gs = zeros(size(g));
                for i = 1:min(numel(g), p.nParam)
                    single_g = g(i);
                    single_pval = pval(i);
                    if i <= numel(p.individualScaling)
                        scaling_type = p.individualScaling{i};
                        box_lims = p.individualBoxLims(i,:);
                    else
                        scaling_type = p.scaling;
                        box_lims = p.boxLims(min(i, size(p.boxLims,1)),:);
                    end
                    gs(i) = p.scaleGradientSingleParameter(single_g, single_pval, scaling_type, box_lims,i);
                end
            elseif ~isempty(p.gradientScalingFunction)
                % Use custom gradient scaling
                gs = p.gradientScalingFunction(g, pval, p);
            else
                % Use built-in gradient scaling
                gs = g.*diff(p.boxLims, [], 2);
                if ~strcmp(p.scaling, 'linear')
                    tmp = (pval-p.boxLims(:,1))./diff(p.boxLims, [], 2);
                    if strcmp(p.scaling, 'exp')
                        gs = gs./dExpScale(tmp, p, []);
                    elseif(strcmp(p.scaling, 'log'))
                        gs = gs./dLogScale(tmp, p, []);
                    end
                end
            end
        end
        
        %------------------------------------------------------------------
        % NEW: Individual scaling methods for single parameters
        function vs = scaleSingleParameter(p, pval, scaling_type, box_lims,i)
            switch scaling_type
                case 'linear'
                    vs = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));
                case 'log'
                    vs = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));
                    vs = logScale(vs,p,i);
                case 'exp'
                    vs = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));
                    vs = expScale(vs,p,i);
                otherwise
                    vs = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));
            end
        end
        
        function pval = unscaleSingleParameter(p, vs, scaling_type, box_lims, i)
            switch scaling_type
                case 'linear'
                    pval = vs * (box_lims(2) - box_lims(1)) + box_lims(1);
                case 'log'
                    vs = expScale(vs, p,i);                                        
                    pval = vs * (box_lims(2) - box_lims(1)) + box_lims(1);
                case 'exp'
                    vs = logScale(vs, p,i);
                    pval = vs * (box_lims(2) - box_lims(1)) + box_lims(1);
                otherwise
                    pval = vs * (box_lims(2) - box_lims(1)) + box_lims(1);
            end
        end
        
        function gs = scaleGradientSingleParameter(p, g, pval, scaling_type, box_lims,i)
            switch scaling_type
                case 'linear'
                    gs = g .* (box_lims(2) - box_lims(1));
                case 'log'              
                    tmp = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));                 
                    gs = g .* (box_lims(2) - box_lims(1));
                    gs = gs ./dLogScale(tmp,p,i);
                case 'exp'
                   tmp = (pval - box_lims(1)) / (box_lims(2) - box_lims(1));
                   gs = g .* (box_lims(2) - box_lims(1));
                    gs = gs ./dExpScale(tmp,p,i);
                otherwise
                    gs = g * (box_lims(2) - box_lims(1));
            end
        end
        
        % ... (rest of your existing methods remain the same)
        %------------------------------------------------------------------
        function g = collapseGradient(p, g)
            assert(strcmp(p.belongsTo, 'state0'),['This function is only',...
                'intended to collapse gradient arrays from parameters that',...
                'belongs to state0'])
            % take sum of each lump
            if ~isempty(p.lumping) && isnumeric(p.lumping)
                g = accumarray(p.lumping, g(p.subset), [], @sum); 
            end
        end
        
        function v = getParameter(p, setup)
            if ~strcmp(p.type, 'multiplier')
                v = p.getParameterValue(setup);
            else
                v = p.getParameterValue(setup, false)./p.referenceValue;
                v = collapseLumps(v, p.lumping);
            end
        end
        
        function setup = setParameter(p, setup, v)
            if ~strcmp(p.type, 'multiplier')
                setup = p.setParameterValue(setup, v);
            else
                v  = expandLumps(v, p.lumping).*p.referenceValue;
                setup = p.setParameterValue(setup, v, false);
            end
        end
        
        function v = getParameterValue(p, setup, doCollapse)
            if nargin < 3
                doCollapse = true;
            end
            if ~(strcmp(p.belongsTo, 'schedule') && strcmp(p.location{1}, 'control'))
                v = p.getfun(setup.(p.belongsTo), p.location{:});
                v = v(p.subset);
                if doCollapse
                    v = collapseLumps(v, p.lumping);
                end
            else % control-parameter (assume constant selection of control steps)
                control = setup.schedule.control(p.controlSteps(1));
                v  = p.getControlParameterValue(control, doCollapse);
            end
        end
        
        function setup = setParameterValue(p, setup, v, doExpand)
            if nargin < 4
                doExpand = true;
            end
            if doExpand
                v  = expandLumps(v, p.lumping);
            end
            if ~(strcmp(p.belongsTo, 'schedule') && strcmp(p.location{1}, 'control'))
                if isnumeric(p.subset)
                    tmp = p.getfun(setup.(p.belongsTo), p.location{:});
                    v   = setSubset(tmp, v, p.subset);
                end
                setup.(p.belongsTo) = ...
                    p.setfunApply(setup.(p.belongsTo), v);
            else % control-parameter (assume constant over selected control steps)
                nc = numel(p.controlSteps);
                for k = 1:nc
                    step = p.controlSteps(k);
                    setup.schedule.control(step) = ...
                        p.setControlParameterValue(setup.schedule.control(step), v);
                end
            end
        end
        %------------------------------------------------------------------       
        function v = getControlParameterValue(p, control, doCollapse)
            if nargin < 3
                doCollapse = true;
            end
            loc = p.location;
            wlocIx = strcmp('W', loc);
            if any(wlocIx)
                wlocIx = find(wlocIx);
                % well-parameter subset applies to well numbers 
                v = applyFunction(@(x)p.getfun(x, loc{wlocIx+1:end}), getfield(control, loc{2:wlocIx}, {p.subset}));
                if ~isempty(p.wellSubsets)
                    v = applyFunction(@(x,y)x(y), v, p.wellSubsets);
                end
                v = vertcat(v{:});
            else
                % apply subset to parameter
                loc = p.location(2:end);
                v = p.getfun(control, loc{:});
                v = v(p.subset);
            end      
            if doCollapse
                v = collapseLumps(v, p.lumping);
            end
        end
        %------------------------------------------------------------------
        function control = setControlParameterValue(p, control, v)
            loc = p.location;
            wlocIx = strcmp('W', loc);
            if any(wlocIx)
                wlocIx = find(wlocIx);
                % well-parameter special treatment
                sub = p.subset;
                W = getfield(control, loc{2:wlocIx});%#ok
                if ~isnumeric(sub)
                    sub = (1:numel(W))';
                end
                % scalar param per well or connection
                if numelValue(v) == numel(sub)
                    np = ones(size(sub));
                else
                    np = arrayfun(@(w)numel(w.cells), W(sub));
                end
                [i1, i2] = deal(cumsum([1;np(1:end-1)]), cumsum(np));
                v  = applyFunction(@(i1,i2)v(i1:i2), i1, i2);
                if ~isempty(p.wellSubsets)
                    for k = 1:numel(sub)
                        tmp  = p.getfun(W(sub(k)), loc{wlocIx+1:end});
                        v{k} = setSubset(tmp, v{k}, p.wellSubset{k});
                    end
                end
                for k = 1:numel(sub)
                    W(sub(k)) = p.setfun(W(sub(k)), loc{wlocIx+1:end}, v{k});
                end
                control = setfield(control, loc{2:wlocIx}, W);%#ok
            else
                % other parameter (set subset of parameter)
                loc = p.location(2:end);
                if isnumeric(p.subset)
                    tmp = p.getfun(control, loc{:});
                    v   = setSubset(tmp, v, p.subset);
                end
                control = p.setfun(control, loc{:}, v);
            end
        end
        %------------------------------------------------------------------
        function s = setfunApply(p, s, v)
            s = p.setfun(s, p.location{:}, v);
            for k = 1:numel(p.extraLocations)
                s = p.setfun(s, p.extraLocations{k}{:}, v);
            end
        end
    end
    
       methods (Access = private)
        function p = setParameterType(p)
            % Set parameter type based on name pattern
            name_lower = lower(p.name);
            if contains(name_lower, {'ne_bg', 'pe_bg', 'elyte_bg'})
                p.parameterType = 'background_conductivity';
            elseif contains(name_lower, {'ne_kappa', 'pe_kappa', 'elyte_kappa'})
                p.parameterType = 'kappa_conductivity';
            elseif contains(name_lower, {'elyte_d', 'ne_d', 'pe_d'})
                p.parameterType = 'diffusivity';
            elseif contains(name_lower, 'equivalent_eff_cond')
                p.parameterType = 'equivalent_conductivity';
            elseif contains(name_lower, 'current_collector')
                p.parameterType = 'current_collector_conductivity';
            else
                p.parameterType = 'generic';
            end
        end
        
        function p = setupCustomScaling(p)
            % Setup custom scaling functions based on parameter type
            switch p.parameterType
                case 'background_conductivity'
                    p.scaling = 'log';
                case 'kappa_conductivity'
                    p.scaling = 'log';
                case 'diffusivity'
                    p.scaling = 'log';
                case 'equivalent_conductivity'
                    p.scaling = 'log';
                case 'current_collector_conductivity'
                    p.scaling = 'linear';
                otherwise
                    p.scaling = 'log';
            end
        end
        
        function p = setupIndividualParameterScaling(p)
            % Setup individual scaling based on parameter shortnames
            % Assumes p.boxLims is always provided and never empty

            if isempty(p.shortnames) || p.nParam ~= numel(p.shortnames)
                p.individualScaling = [];
                p.individualBoxLims = [];
                p.parameterNames = {};
                return;
            end

            p.parameterNames = p.shortnames;

            % Define scaling types based on parameter name patterns
            scalingMap = struct(...
                'bg', 'linear', ...      % Bruggeman coefficients (exponents)
                'vsa', 'log', ...        % Surface areas (wide range)
                'D', 'log', ...          % Diffusion coefficients (orders of magnitude)
                'kappa', 'log', ...      % Conductivities (wide range)
                'elyte', 'log' ...       % Electrolyte parameters
                );

            p.individualScaling = cell(size(p.parameterNames));

            for i = 1:numel(p.parameterNames)
                param_name = p.parameterNames{i};

                if contains(param_name, '_bg')
                    p.individualScaling{i} = scalingMap.bg;
                elseif contains(param_name, '_vsa')
                    p.individualScaling{i} = scalingMap.vsa;
                elseif contains(param_name, '_D')
                    p.individualScaling{i} = scalingMap.D;
                elseif contains(param_name, '_kappa')
                    p.individualScaling{i} = scalingMap.kappa;
                else
                    p.individualScaling{i} = 'log'; % Default
                end
            end

            % Use provided boxLims - handle different sizes
            if size(p.boxLims, 1) == 1
                p.individualBoxLims = repmat(p.boxLims, numel(p.parameterNames), 1);
            else
                p.individualBoxLims = p.boxLims;
            end

            fprintf('Individual scaling enabled for %d parameters:\n', p.nParam);
            for i = 1:numel(p.parameterNames)
                fprintf('  %s: %s scaling\n', p.parameterNames{i}, p.individualScaling{i});
            end
        end
    end
end

%--------------------------------------------------------------------------
% Custom scaling functions for different parameter types
%--------------------------------------------------------------------------
function vs = scaleBackgroundConductivity(pval, p)
    % Logarithmic scaling for background conductivities
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    % Map to [0,1] range using logarithmic scaling
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    vs = (log10(normalized) - log_min) / (log_max - log_min);
end

function pval = unscaleBackgroundConductivity(vs, p)
    % Unscaling for background conductivities
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    normalized = 10.^(vs * (log_max - log_min) + log_min);
    
    if strcmp(p.type, 'multiplier')
        pval = normalized .* p.referenceValue;
    else
        pval = normalized;
    end
end

function gs = scaleGradientBackgroundConductivity(g, pval, p)
    % Gradient scaling for background conductivities
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    gs = g .* normalized * (log_max - log_min) * log(10);
end

function vs = scaleKappaConductivity(pval, p)
    % Custom scaling for kappa parameters - can be adjusted based on your needs
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    % Use logarithmic scaling for kappa as well
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    vs = (log10(normalized) - log_min) / (log_max - log_min);
end

function pval = unscaleKappaConductivity(vs, p)
    % Custom unscaling for kappa parameters
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    normalized = 10.^(vs * (log_max - log_min) + log_min);
    
    if strcmp(p.type, 'multiplier')
        pval = normalized .* p.referenceValue;
    else
        pval = normalized;
    end
end

function gs = scaleGradientKappaConductivity(g, pval, p)
    % Gradient scaling for kappa parameters
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    gs = g .* normalized * (log_max - log_min) * log(10);
end

function vs = scaleEquivalentConductivity(pval, p)
    % Scaling for equivalent effective conductivity
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    % Logarithmic scaling for equivalent conductivity
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    vs = (log10(normalized) - log_min) / (log_max - log_min);
end

function pval = unscaleEquivalentConductivity(vs, p)
    % Unscaling for equivalent effective conductivity
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    normalized = 10.^(vs * (log_max - log_min) + log_min);
    
    if strcmp(p.type, 'multiplier')
        pval = normalized .* p.referenceValue;
    else
        pval = normalized;
    end
end

function gs = scaleGradientEquivalentConductivity(g, pval, p)
    % Gradient scaling for equivalent effective conductivity
    if strcmp(p.type, 'multiplier')
        ref = p.referenceValue;
        normalized = pval ./ ref;
    else
        normalized = pval;
    end
    
    log_min = log10(min(p.boxLims(:,1)));
    log_max = log10(max(p.boxLims(:,2)));
    gs = g .* normalized * (log_max - log_min) * log(10);
end

%--------------------------------------------------------------------------
% Existing helper functions (keep your original implementations)
%--------------------------------------------------------------------------
function p = setupDefaults(p, setup, opt)
% Make sure setup makes sense and add boxLims if not provided
if isempty(p.subset)
    p.subset = ':';
end
if islogical(p.subset)
    p.subset = find(p.subset);
end
for k = 1:numel(p.wellSubsets)
    if islogical(p.wellSubsets{k})
        p.wellSubsets{k} = find(p.wellSubsets{k});
    end
end

if strcmp(p.belongsTo, 'schedule')
    if isempty(p.controlSteps)
        p.controlSteps = (1:numel(setup.schedule.control));
    else
        % set name depending on first selected control step
        p.name = sprintf('%s_step%d', p.name, p.controlSteps(1));
    end
    if iscell(p.lumping)
        p = handleWellLumping(p, setup);
    end
end

v = getParameterValue(p, setup, false);

if ~isempty(p.lumping)
    p.nParam = max(double(p.lumping));
    if numel(p.lumping) == 1 && p.lumping == 1
        p.lumping = ones(numel(v), 1);
    end
else
    p.nParam = numel(v);
end
   
if strcmp(p.type, 'multiplier')
    if any(v==0)
        error('Parameters (''%s'') of type ''multiplier'' can''t be zero.', p.name);
    end
    p.referenceValue = v;
end

% Handle default parameter limits with parameter-specific defaults
if isempty(p.boxLims)
    switch p.parameterType
        case {'background_conductivity', 'kappa_conductivity'}
            p.boxLims = [1e-6, 1e2]; % Appropriate range for conductivities
        case 'diffusivity'
            p.boxLims = [1e-12, 1e-4]; % Typical diffusivity range
        case 'equivalent_conductivity'
            p.boxLims = [1e-3, 1e3]; % Range for equivalent conductivity
        case 'current_collector_conductivity'
            p.boxLims = [1e-6, 1e8]; % Wide range for current collectors
        otherwise
            rlim  = opt.relativeLimits;
            if strcmp(p.type, 'value')
                if opt.uniformLimits
                    tmp = [min(min(v)), max(max(v))];
                    ii  = [1+(tmp(1)<0), 2-(tmp(2)<0)];
                    p.boxLims = tmp.*rlim(ii);
                else
                    p.boxLims = bsxfun(@times, v, rlim);
                end
                isNeg = p.boxLims(:,1) < 0;
                p.boxLims(isNeg,:) = p.boxLims(isNeg, [2, 1]);
            else
                p.boxLims = rlim;
            end
    end
    
    % special treatment of saturations
    if any(strcmp(p.name, {'sw', 'sg'}))
         p.boxLims = [0, 1];
    elseif any(v==0) 
    % check if relative limits are set for zero-value params
        p.boxLims(v==0, 1) = -1; 
        p.boxLims(v==0, 2) =  1;
        warning('Parameter %s contains zero-values. %s',  p.name, ...
         'Defaulting lower/upper limits (''boxLims'') to [-1 1]');
    end
end

assert(any(size(p.boxLims,1) == [1, p.nParam]), ...
    'Property ''boxLims'' does not match number of parameters');
end

%--------------------------------------------------------------------------
function p = setupByName(p, SimulatorSetup)
% setup for typical parameters
[getfun, setfun] = deal([]);
switch lower(p.name)
    case 'transmissibility'
        belongsTo = 'model';
        location = {'operators', 'T'};
    case {'permx', 'permy', 'permz'}
        belongsTo = 'model';
        col = find(strcmpi(p.name(end), {'x', 'y', 'z'}));
        location = {'rock', 'perm', {':', col}};
        setfun   = @setPermeabilityFun;
    case 'porevolume'
        belongsTo = 'model';
        location = {'operators', 'pv'};
    case 'conntrans'
        belongsTo = 'schedule';
        location = {'control', 'W', 'WI'};
    case {'sw', 'sg'}
        belongsTo = 'state0';
        col = SimulatorSetup.model.getPhaseIndex(upper(p.name(end)));
        location = {'s', {':', col}};
        oix = SimulatorSetup.model.getPhaseIndex('O');
        assert(~isempty(oix), ...
            'Current assumption is that oil is the dependent phase');
        setfun   = @(obj, varargin)setSaturationFun(obj, varargin{:}, oix);
    case 'pressure'
        belongsTo = 'state0';
        location = {'pressure'};
    case {'swl', 'swcr', 'swu', 'sowcr', 'sogcr', 'sgl', 'sgcr', ...
            'sgu', 'krw', 'kro', 'krg'}
        belongsTo = 'model';
        map = getScalerMap();
        ix  = map.kw.(upper(p.name));
        if ~iscell(ix), ix = {ix}; end
        [ph, col] = deal(map.ph{ix{1}(1)}, ix{1}(2));
        location = {'rock', 'krscale', 'drainage', ph, {':', col}};
        for k = 2:numel(ix)
            [ph, col] = deal(map.ph{ix{k}(1)}, ix{k}(2));
            p.extraLocations{k-1} = {'rock', 'krscale', 'drainage', ph, {':', col}};
        end
        setfun   = @setRelPermScalersFun;
    case {'bhp', 'rate', 'wrat', 'orat', 'grat'}
        belongsTo = 'schedule';
        location = {'control', 'W', 'val'};
        getfun = @(obj, varargin)getWellControlValue(obj, varargin{:}, ...
                                                     lower(p.name));
        setfun = @(obj, varargin)setWellControlValue(obj, varargin{:}, ...
                                                     lower(p.name));
        p.controlType = lower(p.name);
    otherwise
        error('No default setup for parameter: %s\n', p.name);     
end
if isempty(p.belongsTo)
    p.belongsTo = belongsTo;
end
if isempty(p.location)
    p.location = location;
end
if isempty(p.setfun) && ~isempty(setfun)
    p.setfun = setfun;
end
if isempty(p.getfun) && ~isempty(getfun)
    p.getfun = getfun;
end
end

%--------------------------------------------------------------------------            
function map = getScalerMap()
phOpts = {'w', 'ow', 'g', 'og'};
kw  = struct('SWL',   [1,1], 'SWCR',  [1,2], 'SWU', [1,3], ...
             'SGL',   [3,1], 'SGCR',  [3,2], 'SGU', [3,3], ...
             'SOWCR', [2,2], 'SOGCR', [4,2], ...
             'KRW',   [1,4], 'KRO',   {{[2, 4],[4,4]}}, 'KRG', [3,4]);
map = struct('ph', {phOpts}, 'kw', kw);
end

%--------------------------------------------------------------------------
function v = collapseLumps(v, lumps)
% take mean of each lump
if ~isempty(lumps) && isnumeric(lumps)
    if isa(v, 'double')
        v = accumarray(lumps, v, [], @mean);
    else % special treatment in case of ADI
        M = sparse(lumps, (1:numel(lumps))', 1);
        v = (M*v)./sum(M,2);
    end
end
end

%--------------------------------------------------------------------------
function v = expandLumps(v, lumps)
if ~isempty(lumps) && isnumeric(lumps)
    v = v(lumps);
end
end

%--------------------------------------------------------------------------
function v = setSubset(v, vi, sub)
if isa(vi, 'ADI') && ~isa(v, 'ADI')
    v = double2ADI(v, vi);
end
v(sub) = vi;
end

%--------------------------------------------------------------------------
function checkSetup(p, setup)
if ~isempty(p.lumping)
     tmp = p.getParameterValue(setup, false);
     assert(numel(p.lumping) == numel(tmp), 'Lumping vector has incorrect size');
end
tmp = p.getParameter(setup);
assert(numel(tmp)==p.nParam, 'Report error to develolper')
assert(any(size(p.boxLims,1) == [1, p.nParam]), ...
       'The number of upper/lower limits should be 1 or nParam');
chk = tmp >= p.boxLims(:,1)-1.0e-6 & tmp <= p.boxLims(:,2)+1.0e-6;
assert(all(chk(~isnan(tmp))), 'Parameter values are not within given limits');
if ~strcmp(p.controlType, 'none')
    assert(any(strcmp(p.controlType, {'bhp', 'rate', 'wrat', 'orat', 'grat', 'policy'})), ...
        'Well control of type "%s" is not supported', p.controlType);
end
if ~isempty(p.extraLocations)
    assert(~strcmp(p.belongsTo, 'schedule'), ...
        'Multiple locations for schedule-parameters are not supported');
end
end

%--------------------------------------------------------------------------    
function model = setPermeabilityFun(model, varargin)
% utility for setting permx/y/z possibly as AD and include effect on
% transmissibilities
[location, v] = deal(varargin(1:end-1), varargin{end});
[nc, nd] = size(model.rock.perm);
perm = model.rock.perm;
if ~iscell(perm)
    perm = mat2cell(perm, nc, ones(1,nd));
end
col = location{end}{end};
assert(col<=nd, 'Can''t get column %d since perm has %d column(s)', col, nd);
perm{col} = v;
% transmissibilities
th = 0;
for k = 1:nd
    th = th + perm2directionalTrans(model, perm{k}, k);
end
cf = model.G.cells.faces(:,1);
nf = model.G.faces.num;
% mapping from from cell-face to face
M = sparse(cf, (1:numel(cf))', 1, nf, numel(cf));
% consider only internal faces
ie = model.operators.internalConn;
model.operators.T = 1./(M(ie,:)*(1./th));
% if all perms are doubles, set to nc x nd array
if all(cellfun(@(p)isa(p, 'double'), perm))
    model.rock.perm = cell2mat(perm);
else
    model.rock.perm = perm;
end
end

%--------------------------------------------------------------------------    
function state = setSaturationFun(state, varargin)
assert(nargin >= 3, 'Insufficient input')
[loc, v, oix] = deal(varargin(1:end-2), varargin{end-1}, varargin{end});
assert(isa(v, 'double'), 'Setting saturation to class %s is not supported', class(v));
pix = loc{end}{end};
ds = v-state.s(:, pix);
state.s(:, pix) = v;
state.s(:, oix) =  state.s(:, oix) - ds;
end

%--------------------------------------------------------------------------   
function model = setRelPermScalersFun(model, varargin)
[loc, v] = deal(varargin(1:end-1), varargin{end});
if ~isa(v, 'ADI')
    model = setfield(model, loc{:}, v);
else
    % last location is column no, second last is phase
    col = loc{end}{end};
    ph  = loc{end-1};
    d   = getfield(model, loc{1:end-2});  %#ok
    if ~isfield(d, 'tmp') || ~isfield(d.tmp, ph)
        nc = model.G.cells.num;
        d.tmp.(ph) = mat2cell(d.(ph), nc, ones(1,4));
    end
    d.tmp.(ph){col} = v;
    d.(ph) = @(cells, col)d.tmp.(ph){col}(cells);
    model.rock.krscale.drainage = d;
end
end

%--------------------------------------------------------------------------   
function v = getWellControlValue(w, varargin)
% don't use location
cntr = varargin{end};
if strcmp(w.type, cntr)
    v = w.val;
else
    %return limit-value if exists
    if isfield(w, 'lims') && isfield(w.lims, cntr)
        v = w.lims.(cntr);
    else
        error('Request to set non-existing control/limit %s for well %s.', cntr, w.name);
    end
end
end

%--------------------------------------------------------------------------   
function w = setWellControlValue(w, varargin)
% don't use location
[v, cntr] = deal(varargin{end-1}, varargin{end});
ok = false;
% set both value and limit fields
if strcmp(w.type, cntr)
    w.val = v;
    ok = true;
end
if isfield(w, 'lims') && isfield(w.lims, cntr)
    w.lims.(cntr) = v;
    ok = true;
end
if ~ok
    error('Request to set non-existing control/limit %s for well %s.', cntr, w.name);
end
end

%--------------------------------------------------------------------------          
function ti = perm2directionalTrans(model, p, cdir)
% special utility function for calculating transmissibility along coordinate direction cdir
assert(size(p,2)==1, ...
       'Input p should be single column representing permeability in direction cdir');
% make perm represent diag perm tensor with zero perm orthogonal to cdir
dp = value(p);
r.perm = zeros(numel(dp), 3);
r.perm(:, cdir) = dp;
ti = computeTrans(model.G, r);
if isa(p, 'ADI')
    % make ti ADI (note ti is linear function of p)
    p = p./dp;
    cellno = rldecode(1 : model.G.cells.num, diff(model.G.cells.facePos), 2).';
    ti = ti.*p(cellno);
end
end

%--------------------------------------------------------------------------       
function vs = expScale(v, p, i)
base = getScalingBase(p);
if isempty(i)
    vs = (base.^v - 1)./(base - 1);
else
    vs = (base(i).^v - 1)./(base(i) - 1);
end
end

%--------------------------------------------------------------------------       
function ds = dExpScale(v, p,i)
base = getScalingBase(p);
if isempty(i)
ds = (base.^v).*(log(base)./(base - 1));
else
ds = (base(i).^v).*(log(base(i))./(base(i) - 1));
end
end

%--------------------------------------------------------------------------       
function vs = logScale(v, p,i)
base = getScalingBase(p);
if isempty(i)
    vs = log((base-1).*v+1)./log(base);
else
    vs = log((base(i)-1).*v+1)./log(base(i));
end
end

%--------------------------------------------------------------------------       
function ds = dLogScale(v, p,i)
base = getScalingBase(p);
if isempty(i)
    ds = (base-1)./( ((base-1).*v+1).*log(base) );
else
    ds = (base(i)-1)./( ((base(i)-1).*v+1).*log(base(i)) );
end
end

%--------------------------------------------------------------------------     
function base = getScalingBase(p)
base = p.scalingBase;
if isnan(base)
    base = p.boxLims(:,2)./p.boxLims(:,1);
end
end

%--------------------------------------------------------------------------     
function p = handleWellLumping(p, setup)
% This is a bit intricate
assert(strcmp(p.belongsTo, 'schedule'), 'Expected parameter to be of schedule-type');
W = setup.schedule.control(1).W;
ss = p.subset;
if ~isnumeric(ss)
    ss = (1:numel(W))';
end
if iscell(p.lumping) % lumping given for each well
    assert(numel(ss) == numel(p.lumping), 'Well subset vs well lumping mismatch');
    if any(cellfun(@islogical, p.lumping))
        % a logical (true) is treated as individual well lumping
        assert(all(cellfun(@islogical, p.lumping)));
        p.lumping = num2cell((1:numel(p.lumping)));
    end
    for k = 1:numel(ss)
        nc = numel(W(ss(k)).cells);
        if ~isempty(p.wellSubsets)
            nc = numel(p.wellSubsets{k});
        end
        if numel(p.lumping{k}) == 1 && nc > 1
            p.lumping{k} = ones(nc,1)*p.lumping{k};
        end
    end
    p.lumping = vertcat(p.lumping{:});
end
end

%--------------------------------------------------------------------------
function result = applyFunction(fun, varargin)
% Apply function to multiple inputs
result = fun(varargin{:});
end

%--------------------------------------------------------------------------
function n = numelValue(v)
% Get number of elements in value (handles ADI variables)
if isa(v, 'ADI')
    n = numel(v.val);
else
    n = numel(v);
end
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU
%}