%% Script to plot the results from the high-rate calibration

clear all
close all force

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

getTime = @(states) cellfun(@(s) s.time, states);
getEndTime = @(states) states{end}.time;
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
getI = @(states) cellfun(@(s) s.(ctrl).I, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

tag = 'no-elyte-params';
tag = 'one-elyte-param';
tag = 'two-elyte-params';
% tag = 'three-elyte-params';

diary(sprintf('_diary-%s-%s-%s.txt', mfilename, tag, datestr(now, 'yyyymmdd-HHMMSS')));
dosave = true;


% Load

% Original data
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Highest DRate is last
k = numel(dataraw.time);
expdata = struct('time', dataraw.time{k} * hour, ...
                 'U'   , dataraw.voltage{k}    , ...
                 'I'   , abs(mean(dataraw.current{k})));

Dne = 1e-14;
Dpe = 1e-14;
output14 = load(sprintf('high-rate-calibrated-outputOpt-%s-%g-%g.mat', tag, Dne, Dpe));

Dne = 1e-13;
Dpe = 1e-14;
output13 = load(sprintf('high-rate-calibrated-outputOpt-%s-%g-%g.mat', tag, Dne, Dpe));

%% Convergence reason

r1 = getReasonStr(output14.history);
disp('1e-14:')
disp(r1)
r2 = getReasonStr(output13.history);
disp('1e-13:')
disp(r2)

%% Plot

colors = lines(4);
fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
title(tag);
hold on;
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 2C');

Dne = output14.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
Dpe = output14.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', Dne, Dpe);
plot(getTime(output14.output0.states)/hour, getE(output14.output0.states), 'color', colors(1,:), 'displayname', dn);

Dne = output13.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
Dpe = output13.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Initial guess D_{NE}=%g D_{PE}=%g', Dne, Dpe);
plot(getTime(output13.output0.states)/hour, getE(output13.output0.states), 'color', colors(2,:), 'displayname', dn);

plot(getTime(output14.outputOpt.states)/hour, getE(output14.outputOpt.states), '--', 'color', colors(1,:), 'displayname', 'Calibrated');
plot(getTime(output13.outputOpt.states)/hour, getE(output13.outputOpt.states), '--', 'color', colors(2,:), 'displayname', 'Calibrated');

xlabel('Time  /  h')
ylabel('E  /  V')

legend()
drawnow

%dosave = false;
if dosave
    exportgraphics(fig, sprintf('/tmp/high-rate-calibration-results-%s.png', tag), 'Resolution', 300);
end

%% Plot difference between experiment and calibrated

fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on;
title(tag)

expdataUinterp1 = @(t) interp1(expdata.time, expdata.U, t, 'linear', 'extrap');

% quantify difference
wL214 = sqrt(trapz(getTime(output14.outputOpt.states), (getE(output14.outputOpt.states) - expdataUinterp1(getTime(output14.outputOpt.states))).^2)) / getEndTime(output14.outputOpt.states) / milli;
wL213 = sqrt(trapz(getTime(output13.outputOpt.states), (getE(output13.outputOpt.states) - expdataUinterp1(getTime(output13.outputOpt.states))).^2)) / getEndTime(output13.outputOpt.states) / milli;
fprintf('L2 error Dne=1e-14: %g mV\n', wL214);
fprintf('L2 error Dne=1e-13: %g mV\n', wL213);

Dne = output14.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
Dpe = output14.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Difference with D_{NE}=%g (wL2=%g mV)', Dne, wL214);
plot(getTime(output14.outputOpt.states)/hour, abs(getE(output14.outputOpt.states) - expdataUinterp1(getTime(output14.outputOpt.states))), 'color', colors(1,:), 'displayname', dn);

Dne = output13.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient;
Dpe = output13.output0.model.(pe).(co).(am).(sd).referenceDiffusionCoefficient;
dn = sprintf('Difference with D_{NE}=%g (wL2=%g mV)', Dne, wL213);
plot(getTime(output13.outputOpt.states)/hour, abs(getE(output13.outputOpt.states) - expdataUinterp1(getTime(output13.outputOpt.states))), 'color', colors(2,:), 'displayname', dn);

xlabel('Time  /  h')
ylabel('|E_{exp} - E_{sim}|  /  V')

legend();
drawnow

%dosave = false;
if dosave
    exportgraphics(fig, sprintf('/tmp/high-rate-calibration-results-%s-difference.png', tag), 'Resolution', 300);
end



%% Plot dqdv: may need spline approx of OCVs

doplot = false;
if doplot

    fig = figure; %('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
    hold on;
    expQ = cumtrapz(expdata.time, expdata.I*ones(size(expdata.time))) / hour;

    dqdvExp = gradient(expQ, expdata.U);
    plot(expdata.U, -dqdvExp, 'k-', 'displayname', 'Experiment 2C');

    simQ = @(states) cumtrapz(getTime(states), getI(states)) / hour;

    dqdv = gradient(simQ(output14.outputOpt.states), getE(output14.outputOpt.states));
    plot(getE(output14.outputOpt.states), -dqdv, 'color', colors(1,:), 'displayname', 'Calibrated D_{NE}=1e-14');

    dqdv = gradient(simQ(output13.outputOpt.states), getE(output13.outputOpt.states));
    plot(getE(output13.outputOpt.states), -dqdv, 'color', colors(2,:), 'displayname', 'Calibrated D_{NE}=1e-13');

    xlabel 'E  /  V'
    ylabel '-dQdV  /  Ah V^{-1}'

    legend();

end

%% Plot dvdq: may need spline approx of OCVs

doplot = false;
if doplot

    fig = figure;%('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
    hold on;
    dvdqExp = gradient(expdata.U, expQ);
    plot(expQ, -dvdqExp, 'k--', 'displayname', 'Experiment 2C');

    dvdq = gradient(getE(output14.outputOpt.states), simQ(output14.outputOpt.states));
    plot(simQ(output14.outputOpt.states), -dvdq, 'color', colors(1,:), 'displayname', 'Calibrated D_{NE}=1e-14');
    dvdq = gradient(getE(output13.outputOpt.states), simQ(output13.outputOpt.states));
    plot(simQ(output13.outputOpt.states), -dvdq, 'color', colors(2,:), 'displayname', 'Calibrated D_{NE}=1e-13');
    ylim([0, max(abs(dvdqExp))*0.2]);

    xlabel 'Q  /  Ah'
    ylabel '-dVdQ  /  V A^{-1} h^{-1}'
    legend();

end

%% Calculate tortuosity

tau1init = tortuosity(output14.output0.model);
tau1opt = tortuosity(output14.outputOpt.model);

tau2init = tortuosity(output13.output0.model);
tau2opt = tortuosity(output13.outputOpt.model);

assert(tau1init.(ne) == tau2init.(ne));
assert(tau1init.(pe) == tau2init.(pe));

% Add to jsonstructHRC
output14.jsonstructHRC.(ne).tortuosity = tau1opt.(ne);
output14.jsonstructHRC.(pe).tortuosity = tau1opt.(pe);
output13.jsonstructHRC.(ne).tortuosity = tau2opt.(ne);
output13.jsonstructHRC.(pe).tortuosity = tau2opt.(pe);

% Compare
fjv = jsonDiff(output14.jsonstructHRC, ...
               output13.jsonstructHRC, ...
               sprintf('Dne %g', output14.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient), ...
               sprintf('Dne %g', output13.output0.model.(ne).(co).(am).(sd).referenceDiffusionCoefficient));

% Add column with initial values
model0 = output14.output0.model;
params = fjv.flatjson(:,1);
initvals = cell(size(params));

for iparam = 1:numel(params)
    fn = params{iparam};
    if contains(fn, 'tortuosity')
        if contains(fn, ne)
            initvals{iparam} = tau1init.(ne);
        elseif contains(fn, pe)
            initvals{iparam} = tau1init.(pe);
        else
            error('Unknown electrode in tortuosity param %s', fn);
        end
    else
        fn2 = strsplit(fn, '.');
        initvals{iparam} = getfield(model0, fn2{:});
    end
end

fjv.flatjson = [fjv.flatjson, initvals];
fjv.columnnames = [fjv.columnnames, {'Initial'}];
fjv.flatjson(:,4) = [];
fjv.columnnames(4) = [];
tt = fjv.print();

% shorten names: NegativeElectrode -> NE, PositiveElectrode -> PE
for ir = 1:size(tt,1)
    tt{ir,1} = strrep(tt{ir,1}, 'NegativeElectrode', 'NE');
    tt{ir,1} = strrep(tt{ir,1}, 'PositiveElectrode', 'PE');
    tt{ir,1} = strrep(tt{ir,1}, 'Separator', 'Sep');
    tt{ir,1} = strrep(tt{ir,1}, 'ActiveMaterial', 'AM');
    %tt{ir,1} = strrep(tt{ir,1}, 'SolidDiffusion', 'SD');
    tt{ir,1} = strrep(tt{ir,1}, 'referenceDiffusionCoefficient', 'D');
    tt{ir,1} = strrep(tt{ir,1}, 'volumetricSurfaceArea', 'VSA');
end

fig = figure('Units', 'inches', 'Position', [0.1, 0.1, 11, 4]);
title(tag)
% Get the table in string form.
tstr = evalc('disp(tt)');
tstr = strrep(tstr, '{', '');
tstr = strrep(tstr, '}', '');
tstr = strrep(tstr, '<strong>', '');
tstr = strrep(tstr, '</strong>', '');
fixedwidthfont = get(0, 'FixedWidthFontName');
annotation(gcf, 'Textbox', 'String', tstr, 'Position', [0 -0.1 1 1], ...
           'fontname', fixedwidthfont, ...
           'fontsize', 14, 'interpreter', 'laTex');
axis off
drawnow

if dosave
    exportgraphics(fig, sprintf('/tmp/high-rate-calibration-results-%s-table.png', tag), 'Resolution', 300);
end


% % Use TeX Markup for bold formatting and underscores.
% %tstr = strrep(tstr,'<strong>','\bf');
% %tstr = strrep(tstr,'</strong>','\rm');
% %tstr = strrep(tstr,'_','\_');
% %tstr = strrep(tstr,' ', '\_');
% % latex
% tstr = strrep(tstr, '{', '');
% tstr = strrep(tstr, '}', '');
% tstr = strrep(tstr, '<strong>', '');
% tstr = strrep(tstr, '</strong>', '');
% % Get a fixed-width font.
% FixedWidth = get(0,'FixedWidthFontName');
% % Output the table using the annotation command.
% % annotation(gcf,'Textbox','String',tstr,'Interpreter','Tex',...
% %            'FontName',FixedWidth,'fontsize', 12), 'Units','Normalized','Position',[0 0 1 1]);
% annotation(gcf, 'Textbox', 'String', tstr, 'Position', [0 -0.1 1 1], 'fontname', FixedWidth, ...
%           'fontsize', 12, 'interpreter', 'laTex');
% axis off

% fig = uifigure('Name', tag, 'Position', [100, 100, 800, 400]);
% uit = uitable(fig, "Data", tt);
% uit.Position(3) = fig.Position(3) - 80;
% uit.Position(4) = fig.Position(4) - 80;

drawnow

diary off;
