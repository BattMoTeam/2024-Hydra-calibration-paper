%% Script to calibrate parameters under equilibrium assumptions

clear all
close all

diary(sprintf('diary-%s.txt', datestr(now, 'yyyymmdd-HHMMSS')));

mrstDebug(0);

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 15);
set(0, 'defaultaxesfontsize', 15);

am   = 'ActiveMaterial';
itf  = 'Interface';
pe   = 'PositiveElectrode';
ne   = 'NegativeElectrode';
co   = 'Coating';
sd   = 'SolidDiffusion';
ctrl = 'Control';

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

%% Fetch experimental data

% Smooth voltages
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473smooth');
saveddata    = load(datafilename);
datasmooth   = saveddata.expsmooth;

% Original data (for currents)
datafilename = fullfile(getHydra0Dir(), 'rawData', 'TE_1473.mat');
saveddata    = load(datafilename);
dataraw      = saveddata.experiment;

% Pack data with lowest rate (first item)
expdata = struct('time' , datasmooth.time{1} * hour                           , ...
                 'U'    , datasmooth.voltage{1}                               , ...
                 'cap'  , abs(trapz(dataraw.time{1}*hour, dataraw.current{1})), ...
                 'I'    , abs(mean(dataraw.current{1}))                       , ...
                 'DRate', datasmooth.crate{1});

%% Initial guess simulation

input = struct('DRate'    , expdata.DRate, ...
               'totalTime', expdata.time(end));
outputInit = runHydra(input, 'clearSimulation', false);
css0 = CellSpecificationSummary(outputInit.model);

%% Setup and run optimization

ecs = EquilibriumCalibrationSetup2222(outputInit.model, expdata);
ecs = ecs.setupCalibrationCase(1, 'np_ratio', css0.NPratio);

doipopt = false;

if doipopt
    ipoptOptions = struct('print_level', 5, ...
                          'tol', 1e-5);
    [Xopt, info] = ecs.runIpOpt(ipoptOptions);
    iter = info.iter;
else
    [Xopt, history] = ecs.runUnitBoxBFGS('plotEvolution', false);
    iter = numel(history.val);
end

X0 = ecs.X0;
v0 = ecs.objective(X0);
vopt = ecs.objective(Xopt);
[fexp, fcomp] = ecs.setupfunction();

%% Print

fprintf('\t\tInitial \t Optimized\n');
elde = {'PE', 'NE'};
for e = 1:2
    fprintf('%s\n', elde{e});
    fprintf('theta100\t %1.5f \t %1.5f\n', X0(2*e-1), Xopt(2*e-1));
    fprintf('alpha    \t %1.5f \t %1.5f\n', X0(2*e), Xopt(2*e));
end

fprintf('obj val=%1.2f (%1.2f), iter=%d\n', vopt, v0, iter);

%% Extract parameters

%jsonstructEC = EC.extractAlpha(outputInit.model, Xopt);
ecs.totalAmountVariableChoice = 'volumeFraction';
jsonstructEC = ecs.exportParameters(Xopt);
filename = fullfile(getHydra0Dir(), 'parameters', 'equilibrium-calibration-parameters.json');
writeJsonStruct(jsonstructEC, filename);
printer(jsonstructEC);

%% Simulation with calibrated parameters

input = struct('DRate'        , expdata.I * hour / expdata.cap, ...
               'totalTime'    , expdata.time(end)             , ...
               'lowRateParams', jsonstructEC);
outputOpt = runHydra(input, 'clearSimulation', false);
cssOpt = CellSpecificationSummary(outputOpt.model);
fprintf('NPratio after calibration: %g\n', cssOpt.NPratio);

%% Plot

colors = lines(4);
fig = figure%;('Units', 'inches', 'Position', [0.1, 0.1, 8, 6]);
hold on
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(expdata.time/hour, fcomp(expdata.time, X0), 'color', colors(3,:), 'displayname', 'Initial data');
plot(expdata.time/hour, fcomp(expdata.time, Xopt), 'color', colors(4,:), 'displayname', 'After cell balancing');
%plot(getTime(outputInit.states)/hour, getE(outputInit.states), 'color', colors(1,:), 'displayname', 'P2D initial guess');
plot(getTime(outputOpt.states)/hour, getE(outputOpt.states), 'color', colors(2,:), 'displayname', 'P2D after cell balancing');
xlabel 'Time  /  h';
ylabel 'E  /  V';
legend('location', 'sw')
axis tight
ylim([3.45, 4.9])

dosave = false;
if dosave
    exportgraphics(fig, 'cell-balancing.png', 'resolution', 300);
end

return


%% Plot cell balancing time domain

t = expdata.time;

doextendtime = true;

if doextendtime
    % extend time for the initial before and after
    T = t(end) - t(1);
    ta = t(1); %-0.1*T;
    tb = t(end)+0.4*T;
    t0 = linspace(ta, tb, numel(t));
else
    t0 = t;
end

[~, fpe0, fne0, thetape0, thetane0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne, thetape, thetane] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

if doextendtime
    % truncate PE to cutoff voltage
    cutoff = outputOpt.jsonstruct.(ctrl).lowerCutoffVoltage;
    idx = find(fpe0 <= cutoff, 1, 'first');
    tpe0cut = t0(1:idx);
    fpe0cut = fpe0(1:idx);
    ocp0cut = ocp0(1:idx);

    % truncate NE to max opt voltage
    cutoff = max(fne);
    idx = find(fne0 >= cutoff, 1, 'first');
    tne0cut = t0(1:idx);
    fne0cut = fne0(1:idx);
else
    tpe0cut = t0;
    fpe0cut = fpe0;
    tne0cut = t0;
    fne0cut = fne0;
    ocp0cut = ocp0;
end

figure; hold on; grid on
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(tpe0cut/hour, fpe0cut, 'displayname', 'pe init', 'color', colors(1,:), 'linestyle', '--');
plot(tne0cut/hour, fne0cut, 'displayname', 'ne init', 'color', colors(2,:), 'linestyle', '--');
plot(tpe0cut/hour, ocp0cut, 'displayname', 'ocp init', 'color', colors(3,:), 'linestyle', '--');
plot(t/hour, fpe, 'displayname', 'pe opt', 'color', colors(1,:));
plot(t/hour, fne, 'displayname', 'ne opt', 'color', colors(2,:));
plot(t/hour, ocp, 'displayname', 'ocp opt', 'color', colors(3,:));

xlabel 'Time  /  h';
ylabel 'Voltage  /  V';
legend('location', 'sw')

axis tight
breakyaxis([1, 3]);

%% Plot the same in capacity domain (we know from below that this is
%% just a scaling)

t = expdata.time;
I = expdata.I;

doextendtime = true;
capoverarea = true;
plotdvdq = true;

if doextendtime
    % extend time for the initial before and after
    T = t(end) - t(1);
    ta = t(1); %-0.1*T;
    tb = t(end)+0.4*T;
    t0 = linspace(ta, tb, numel(t))';
else
    t0 = t;
end

q0 = cumtrapz(t0, I*ones(size(t)));
q = cumtrapz(t, I*ones(size(t)));

if capoverarea
    area = outputOpt.jsonstruct.Geometry.faceArea / centi^2; % cm^2
    q0 = q0 / area / milli;
    q = q / area / milli;
    xlab = 'Capacity / mAh{\cdot}cm^{-2}';
else
    xlab = 'Capacity  /  Ah';
end

[~, fpe0, fne0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

if doextendtime
    % truncate PE to cutoff voltage
    cutoff = outputOpt.jsonstruct.(ctrl).lowerCutoffVoltage;
    idx = find(fpe0 <= cutoff, 1, 'first');
    qpe0cut = q0(1:idx);
    fpe0cut = fpe0(1:idx);
    ocp0cut = ocp0(1:idx);

    % truncate NE to max opt voltage
    cutoff = max(fne);
    idx = find(fne0 >= cutoff, 1, 'first');
    qne0cut = q0(1:idx);
    fne0cut = fne0(1:idx);
else
    qpe0cut = q0;
    fpe0cut = fpe0;
    qne0cut = q0;
    fne0cut = fne0;
    ocp0cut = ocp0;
end

figure; hold on; grid on
plot(q/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(qpe0cut/hour, fpe0cut, 'displayname', 'PE init', 'color', colors(1,:), 'linestyle', '--');
plot(qne0cut/hour, fne0cut, 'displayname', 'NE init', 'color', colors(2,:), 'linestyle', '--');
plot(qpe0cut/hour, ocp0cut, 'displayname', 'Cell init', 'color', colors(3,:), 'linestyle', '--');
plot(q/hour, fpe, 'displayname', 'PE opt', 'color', colors(1,:));
plot(q/hour, fne, 'displayname', 'NE opt', 'color', colors(2,:));
plot(q/hour, ocp, 'displayname', 'Cell opt', 'color', colors(3,:));

xlabel(xlab);
ylabel 'Voltage  /  V';
legend('location', 'sw')

axis tight
breakyaxis([1, 3]);

%%

if plotdvdq

    qnom = expdata.cap;
    dvdq = @(v, q) -qnom * gradient(v(:)) ./ gradient(q(:));

    figure; hold on; grid on; legend
    plot(q/hour, dvdq(expdata.U, q0), 'k--', 'displayname', 'Experiment 0.05 C');
    plot(qpe0cut/hour, dvdq(fpe0cut, qpe0cut), 'displayname', 'PE init', 'color', colors(1,:), 'linestyle', '--');
    plot(qne0cut/hour, dvdq(fne0cut, qne0cut), 'displayname', 'NE init', 'color', colors(2,:), 'linestyle', '--');
    plot(qpe0cut/hour, dvdq(ocp0cut, qpe0cut), 'displayname', 'Cell init', 'color', colors(3,:), 'linestyle', '--');
    plot(q/hour, dvdq(fpe, q), 'displayname', 'PE opt', 'color', colors(1,:));
    plot(q/hour, dvdq(fne, q), 'displayname', 'NE opt', 'color', colors(2,:));
    plot(q/hour, dvdq(ocp, q), 'displayname', 'Cell opt', 'color', colors(3,:));
    xlabel(xlab);
    ylabel '-Q_0{\cdot}dV/dQ  /  V';

    ylim([-10, 10])

end


%% Plot cell balancing over capacity (this is just a scaling)

valsInit = ecs.getDefaultValue(); %% !! not ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
[valsOpt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

jsonInit = outputInit.jsonstruct;
jsonOpt = outputOpt.jsonstruct;

minmax = @(x) [min(x), max(x), x(1), x(end)];

% expdata
t = expdata.time;
I = expdata.I;
qexp = cumtrapz(t, I*ones(size(t)));

% pe
mpeInit = valsInit.(pe).totalAmount*ecs.F;
xpeInit = valsInit.(pe).guestStoichiometry100 + qexp / mpeInit;
ypeInit = computeOCPcathodeH0b(xpeInit, [], 1);

mpeOpt = valsOpt.(pe).totalAmount*ecs.F;
xpeOpt = valsOpt.(pe).guestStoichiometry100 + qexp / mpeOpt;
ypeOpt = computeOCPcathodeH0b(xpeOpt, [], 1);

% ne
mneInit = valsInit.(ne).totalAmount*ecs.F;
xneInit = valsInit.(ne).guestStoichiometry100 - qexp / mneInit;
yneInit = computeOCPanodeH0b(xneInit, [], 1);

mneOpt = valsOpt.(ne).totalAmount*ecs.F;
xneOpt = valsOpt.(ne).guestStoichiometry100 - qexp / mneOpt;
yneOpt = computeOCPanodeH0b(xneOpt, [], 1);

% ocp
ocpInit = ypeInit - yneInit;
ocpOpt = ypeOpt - yneOpt;

% plot
colors = lines(4);
figure; hold on; grid on; legend
plot(qexp/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(qexp/hour, ypeInit, '--', 'color', colors(2,:), 'displayname', 'PE init');
plot(qexp/hour, yneInit, '--', 'color', colors(1,:), 'displayname', 'NE init');
plot(qexp/hour, ocpInit, '--', 'color', colors(3,:), 'displayname', 'OCP init');

xlabel 'Capacity  /  Ah';
ylabel 'Voltage  /  V';

%
plot(qexp/hour, ypeOpt, '-', 'color', colors(2,:), 'displayname', 'PE opt');
plot(qexp/hour, yneOpt, '-', 'color', colors(1,:), 'displayname', 'NE opt');
plot(qexp/hour, ocpOpt, '-', 'color', colors(3,:), 'displayname', 'OCP opt');

axis tight
breakyaxis([1, 3]);


%% Plot cell balancing vs capacity calculated using guestStoichiometries: cannot subtract

vals0 = ecs.getDefaultValue(); %% !! not ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
[valsOpt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

jsonInit = outputInit.jsonstruct;
jsonOpt = outputOpt.jsonstruct;

minmax = @(x) [min(x), max(x), x(1), x(end)];

% expdata
t = expdata.time;
I = expdata.I;
qexp = cumtrapz(t, I*ones(size(t)));

% ne
% xne = linspace(vals0.(ne).guestStoichiometry100, vals0.(ne).guestStoichiometry0, 1000)';
% yne = computeOCPanodeH0b(xne, [], 1);

% mne = vals0.(ne).totalAmount*ecs.F;
% xne0 = vals0.(ne).guestStoichiometry100 - qexp / mne;
% disp(minmax(xne0));
% yne0 = computeOCPanodeH0b(xne0, [], 1);
% qne0 = (mne * (vals0.(ne).guestStoichiometry100 - xne0));

mne = vals0.(ne).totalAmount*ecs.F;
xne0 = linspace(vals0.(ne).guestStoichiometry100, vals0.(ne).guestStoichiometry0, 1000)';
disp(minmax(xne0));
yne0 = computeOCPanodeH0b(xne0, [], 1);
qne0 = xne0*mne;
disp(minmax(qne0));

% pe
mpe = vals0.(pe).totalAmount*ecs.F;
%xpe0 = vals0.(pe).guestStoichiometry100 + qexp / mpe;
xpe0 = linspace(vals0.(pe).guestStoichiometry0, vals0.(pe).guestStoichiometry100, 1000)';
disp(minmax(xpe0));
ype0 = computeOCPcathodeH0b(xpe0, [], 1);
%qpe0 = (mpe * (xpe0 - vals0.(pe).guestStoichiometry100));
qpe0 = xpe0*mpe;
disp(minmax(qpe0));

% ocp
ocp = ype0 - yne0;

% plot
colors = lines(4);

figure; hold on; grid on; legend
plot(qexp/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(qne0/hour, yne0, '--', 'color', colors(1,:), 'displayname', 'NE init');
plot(qpe0/hour, ype0, '--', 'color', colors(2,:), 'displayname', 'PE init');
%plot(qocp/hour, ype0 - yne0, '--', 'color', colors(3,:), 'displayname', 'OCP init');

xlabel 'Capacity  /  Ah';
ylabel 'Voltage  /  V';

return

% mne = valsOpt.(ne).totalAmount*ecs.F;
% %xne0 = q / mne + valsOpt.(ne).guestStoichiometry0;
% xneOpt = valsOpt.(ne).guestStoichiometry100 - qexp / mne;
% yneOpt = computeOCPanodeH0b(xneOpt, [], 1);

% mpe = valsOpt.(pe).totalAmount*ecs.F;
% xpeOpt = valsOpt.(pe).guestStoichiometry100 + qexp / mpe;
% ypeOpt = computeOCPcathodeH0b(xpeOpt, [], 1);

% fprintf('opt np ratio: %g\n', mpe/mne);

% colors = lines(4);

% plot(q/hour, yneOpt, '-', 'color', colors(1,:), 'displayname', 'NE Opt');
% plot(q/hour, ypeOpt, '-', 'color', colors(2,:), 'displayname', 'PE Opt');
% plot(q/hour, ypeOpt - yneOpt, '-', 'color', colors(3,:), 'displayname', 'OCP Opt');




% return



%% plot cell balancing capacity domain: only scale x axis by current

t = expdata.time;

% extend time 20% before and after
factor = 0.1;
T = t(end) - t(1);
ta = t(1)-factor*T;
tb = t(end)+factor*T;
t0 = linspace(ta, tb, numel(t));

t0 = t;

% capacities: raw is a straight line (adjust sign)
capraw = -cumtrapz(dataraw.time{1}*hour, dataraw.current{1}); % straight, adjust sign

% processed: plot(cap) is not straight, but plot(expdata.time, cap) is straight
capexp = cumtrapz(expdata.time, ecs.expI*ones(size(expdata.time)));

cap0 = cumtrapz(t0, ecs.expI*ones(size(t0)));
cap = cumtrapz(t, ecs.expI*ones(size(t)));

[~, fpe0, fne0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne, thetape, thetane] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

figure; hold on; grid on
plot(capexp/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(cap0/hour, fpe0, 'displayname', 'pe init', 'color', colors(1,:), 'linestyle', '--');
plot(cap0/hour, fne0, 'displayname', 'ne init', 'color', colors(2,:), 'linestyle', '--');
plot(cap0/hour, ocp0, 'displayname', 'ocp init', 'color', colors(3,:), 'linestyle', '--');
plot(cap/hour, fpe, 'displayname', 'pe opt', 'color', colors(1,:));
plot(cap/hour, fne, 'displayname', 'ne opt', 'color', colors(2,:));
plot(cap/hour, ocp, 'displayname', 'ocp opt', 'color', colors(3,:));

xlabel 'Capacity  /  Ah';
ylabel 'Voltage  /  V';
legend('location', 'sw')

title('Cell balancing in capacity domain');

% % draw vertical lines at exptime limits
% yL = ylim;
% plot(expdata.time(1)*[1, 1]/hour, yL, 'k:', 'handlevisibility', 'off');
% plot(expdata.time(end)*[1, 1]/hour, yL, 'k:', 'handlevisibility', 'off');

%% Cell balancing vs lithiation one electrode at a time to understand better

% Discharge => PE x increases, SOC of battery decrease, OCP decreases

% Plot step by step: first initial, then add optimized

colors = lines(5);
fsz = 14;
style = {'linestyle', ':', 'handlevisibility', 'off'};
left = {'horizontalalignment', 'left', 'fontsize', fsz, 'handlevisibility', 'off'};
right = {'horizontalalignment', 'right', 'fontsize', fsz, 'handlevisibility', 'off'};

%vals0 = ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
%[valsopt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

jsonInit = outputInit.jsonstruct;
jsonOpt = outputOpt.jsonstruct;

t = expdata.time;

% Possibly extend time 20% before and after
factor = 0.1;
T = t(end) - t(1);
ta = t(1)-factor*T;
tb = t(end)+factor*T;
t0 = linspace(ta, tb, numel(t));

t0 = t;

[~, fpeInit, ~, thetaPEinit] = ecs.computeF(t0, X0);
[~, fpeOpt, ~, thetaPEopt, ~] = ecs.computeF(t, Xopt);

figure; hold on; grid on

% 0. Plot original OCPs from xlsx
guestStoichiometry0 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry100;
x = linspace(guestStoichiometry0, guestStoichiometry100, 100);
y = computeOCPcathodeH0b(x, [], 1);
dn = sprintf('PE computeOCP over orig guestStoichiometry (%g, %g)', min(x), max(x));
plot(x, y, 'displayname', dn, 'color', 'k', 'linestyle', ':', 'linewidth', 1)

% include guest stoichiometries start and end
guestStoichiometry0 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry100;
plotline = @(x) line([x, x], [2.9, 5], 'color', 'k', style{:});
plotline(guestStoichiometry0);
plotline(guestStoichiometry100);
plottext = @(x, s) text(x, 4.1, sprintf('PE %s init=%g', s, x), 'color', 'k', left{:});
plottext(guestStoichiometry0, 'guestStoichiometry0');
plottext(guestStoichiometry100, 'guestStoichiometry100');
% text(guestStoichiometry0, 4.1, sprintf('PE guestStoichiometry0 init=%g', guestStoichiometry0), 'color', 'k', left{:});
% text(guestStoichiometry100, 4.1, sprintf('PE guestStoichiometry100 init=%g', guestStoichiometry100), 'color', 'k', left{:});

% 1. Plot initial OCPs in ecs class
color = colors(1,:);
dn = sprintf('fpeinit over thetapeinit (%g, %g)', thetaPEinit(1), thetaPEinit(end));
plot(thetaPEinit, fpeInit, 'color', color, 'linestyle', '--', 'displayname', dn);

% include start and end points
plot(thetaPEinit(1), fpeInit(1), 'o', 'color', color, 'handlevisibility', 'off');
plot(thetaPEinit(end), fpeInit(end), 'x', 'color', color, 'handlevisibility', 'off');
plotline(thetaPEinit(1));
plotline(thetaPEinit(end));
% text(thetaPEinit(1), 3.5, sprintf('PE thetaPEinit(1)=%g', thetaPEinit(1)), 'color', color, left{:});
% text(thetaPEinit(end), 3.5, sprintf('PE thetaPEinit(end)=%g', thetaPEinit(end)), 'color', color, right{:});

fprintf('PE guestStoichiometry0 and guestStoichiometry100 from JsonInit: %1.5f\t %1.5f\n', guestStoichiometry0, guestStoichiometry100);
fprintf('PE guestStoichiometry0 and guestStoichiometry100 from computeF: %1.5f\t %1.5f\n', thetaPEinit(end), thetaPEinit(1));

xlim([-0.1, 1.1]);
ylim([2.9, 5]);

xlabel 'Stoichiometry  /  -';
title('PE OCP vs stoichiometry');
legend('location', 'sw')

return

%% 2a. Plot calibrated OCPs (PE)

% is it like this? or like 2b?

color = colors(2,:);

dn = sprintf('fpeOpt over thetaPEopt (%g, %g)', thetaPEopt(1), thetaPEopt(end));
plot(thetaPEopt, fpeOpt, 'displayname', dn, 'color', color, 'linestyle', ':');

% include start and end points
plot(thetaPEopt(1), fpeOpt(1), 'o', 'color', color, 'handlevisibility', 'off');
plot(thetaPEopt(end), fpeOpt(end), 'x', 'color', color, 'handlevisibility', 'off');

% include guest stoichiometries start and end
guestStoichiometry0 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry100;
dn = sprintf('PE guestStoichiometry0 opt=%g', guestStoichiometry0);
plot(guestStoichiometry0, fpeOpt(1), '.', 'markersize', 12, 'color', color, 'displayname', dn);
plotline = @(x) line([x, x], [4, 5], 'color', color, style{:}, 'handlevisibility', 'off');
plotline(guestStoichiometry0);
plotline(guestStoichiometry100);
text(guestStoichiometry0, 4.5, sprintf('PE guestStoichiometry0 opt=%g', guestStoichiometry0), 'color', color, right{:});
text(guestStoichiometry100, 4.5, sprintf('PE guestStoichiometry100 opt=%g', guestStoichiometry100), 'color', color, right{:});

fprintf('At optimum\n');
fprintf('guestStoichiometry0: %1.5f\n', guestStoichiometry0);
fprintf('guestStoichiometry100: %1.5f\n', guestStoichiometry100);
fprintf('thetaPEopt(1): %1.5f\n', thetaPEopt(1));
fprintf('thetaPEopt(end): %1.5f\n', thetaPEopt(end));

%% 2b. Plot calibrated OCPs with stretch (PE)

guestStoichiometry0 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry100;
x = linspace(guestStoichiometry0, guestStoichiometry100, 1000);
y = computeOCPcathodeH0b(x, [], 1);

guestStoichiometry0 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry100;
xopt = linspace(guestStoichiometry0, guestStoichiometry100, numel(y));

color = colors(2,:);
plot(xopt, y, 'color', color, 'displayname', 'PE computeOCP over opt guestStoichiometry')
plotline = @(x) line([x, x], [4, 5], 'color', color, style{:}, 'handlevisibility', 'off');
plotline(guestStoichiometry0);
plotline(guestStoichiometry100);
text(guestStoichiometry0, 4.5, sprintf('PE guestStoichiometry0 opt=%g', guestStoichiometry0), 'color', color, right{:});
text(guestStoichiometry100, 4.5, sprintf('PE guestStoichiometry100 opt=%g', guestStoichiometry100), 'color', color, right{:});











%% Same for NE (not finished)

% Discharge: NE x decreases, voltage decreases, SOC of battery decreases

colors = lines(5);
fsz = 14;
style = {'linestyle', ':', 'handlevisibility', 'off'};
left = {'horizontalalignment', 'left', 'fontsize', fsz, 'handlevisibility', 'off'};
right = {'horizontalalignment', 'right', 'fontsize', fsz, 'handlevisibility', 'off'};

%vals0 = ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
%[valsopt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

jsonInit = outputInit.jsonstruct;
jsonOpt = outputOpt.jsonstruct;

t = expdata.time;

% Possibly extend time 20% before and after
factor = 0.1;
T = t(end) - t(1);
ta = t(1)-factor*T;
tb = t(end)+factor*T;
t0 = linspace(ta, tb, numel(t));

t0 = t;

[~, ~, fneInit, ~, thetaNEinit] = ecs.computeF(t0, X0);
[~, ~, fneOpt, ~, thetaNEopt] = ecs.computeF(t, Xopt);

figure; hold on; grid on

% 1. Plot initial OCPs
plot(thetaNEinit, fneInit, 'displayname', 'NE init', 'color', colors(3,:), 'linestyle', '--');

% include start and end points
plot(thetaNEinit(1), fneInit(1), 'o', 'color', colors(3,:), 'handlevisibility', 'off');
plot(thetaNEinit(end), fneInit(end), 'x', 'color', colors(3,:), 'handlevisibility', 'off');

% include guest stoichiometries start and end
guestStoichiometry0 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry100;
line([guestStoichiometry0, guestStoichiometry0], [0, 1], 'color', colors(3,:), style{:});
line([guestStoichiometry100, guestStoichiometry100], [0, 1], 'color', colors(3,:), style{:});
text(guestStoichiometry0, 0.75, 'NE guestStoichiometry0 init', 'color', colors(3,:), left{:});
text(guestStoichiometry100, 0.75, 'NE guestStoichiometry100 init', 'color', colors(3,:), left{:});

xlim([-0.3, 1.1]);
ylim([0, 1]);

xlabel 'Stoichiometry \theta  /  -';
title('NE cell balancing vs "lithiation"')
legend;

return

%% 2. Plot calibrated OCPs (NE)

[~, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

plot(thetaNEopt, fneOpt, 'displayname', 'NE opt', 'color', colors(4,:), 'linestyle', ':');

% include start and end points
plot(thetaNEopt(1), fneOpt(1), 'o', 'color', colors(4,:), 'handlevisibility', 'off');
plot(thetaNEopt(end), fneOpt(end), 'x', 'color', colors(4,:), 'handlevisibility', 'off');

% include guest stoichiometries start and end
guestStoichiometry0 = valsNotTruncated.(ne).guestStoichiometry0;
%guestStoichiometry0 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry0;
guestStoichiometry100 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry100;
line([guestStoichiometry0, guestStoichiometry0], [0, 1], 'color', colors(4,:), style{:});
line([guestStoichiometry100, guestStoichiometry100], [0, 1], 'color', colors(4,:), style{:});
text(guestStoichiometry0, 0.75, 'NE guestStoichiometry0 opt', 'color', colors(4,:), right{:});
text(guestStoichiometry100, 0.75, 'NE guestStoichiometry100 opt', 'color', colors(4,:), right{:});



return



colors = lines(5);

vals0 = ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
[valsopt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

t = expdata.time;

% extend time 20% before and after
factor = 0.1;
T = t(end) - t(1);
ta = t(1)-factor*T;
tb = t(end)+factor*T;
t0 = linspace(ta, tb, numel(t));

t0 = t;

[~, fpe0, fne0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne, thetape, thetane] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

figure; hold on; grid on
plot(thetane0, fne0, 'displayname', 'NE init', 'color', colors(2,:), 'linestyle', '--');
plot(thetane, fne, 'displayname', 'NE opt', 'color', colors(4,:), 'linestyle', ':');

% include start and end points
plot(thetane0(1), fne0(1), 'o', 'color', colors(2,:), 'handlevisibility', 'off');
plot(thetane0(end), fne0(end), 'x', 'color', colors(2,:), 'handlevisibility', 'off');
plot(thetane(1), fne(1), 'o', 'color', colors(4,:), 'handlevisibility', 'off');
plot(thetane(end), fne(end), 'x', 'color', colors(4,:), 'handlevisibility', 'off');

% include guest stoichiometries at start and end
plotGuestStoichiometries = true;
if plotGuestStoichiometries
    jsonInit = outputInit.jsonstruct;
    jsonOpt = outputOpt.jsonstruct;

    fsz = 14;
    style = {'linestyle', ':', 'handlevisibility', 'off'};
    left = {'horizontalalignment', 'left', 'fontsize', fsz, 'handlevisibility', 'off'};
    right = {'horizontalalignment', 'right', 'fontsize', fsz, 'handlevisibility', 'off'};

    theta0 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [0, 1], 'color', colors(2,:), style{:});
    line([theta100, theta100], [0, 1], 'color', colors(2,:), style{:});
    text(theta0, 0.1, 'NE theta0 init', 'color', colors(2,:), right{:});
    text(theta100, 0.1, 'NE theta100 init', 'color', colors(2,:), left{:});

    theta0 = valsNotTruncated.(ne).guestStoichiometry0;
    %theta0 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [0, 1], 'color', colors(4,:), style{:});
    line([theta100, theta100], [0, 1], 'color', colors(4,:), style{:});
    text(theta0, 0.75, 'NE theta0 opt', 'color', colors(4,:), right{:});
    text(theta100, 0.75, 'NE theta100 opt', 'color', colors(4,:), left{:});

end

xlabel 'Stoichiometry \theta  /  -';
title('NE cell balancing vs "lithiation"')
legend;



%% Plot cell balancing vs lithiation

% each time gives thetape and thetane, from which we calculate
% ocps. so for each time, we get (thetape, ocppe) and (thetane, ocpne)

% Don't do this plotting, do one electrode at a time

colors = lines(5);

vals0 = ecs.updateGuestStoichiometries(X0, 'includeGuestStoichiometry0', true);
[valsopt, valsNotTruncated] = ecs.updateGuestStoichiometries(Xopt, 'includeGuestStoichiometry0', true);

t = expdata.time;

% extend time 20% before and after
factor = 0.1;
T = t(end) - t(1);
ta = t(1)-factor*T;
tb = t(end)+factor*T;
t0 = linspace(ta, tb, numel(t));

t0 = t;

[~, fpe0, fne0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne, thetape, thetane] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

figure; hold on; grid on
plot(thetape0, fpe0, 'displayname', 'PE init', 'color', colors(1,:), 'linestyle', '--');
plot(thetane0, fne0, 'displayname', 'NE init', 'color', colors(2,:), 'linestyle', '--');
plot(thetape, fpe, 'displayname', 'PE opt', 'color', colors(3,:), 'linestyle', ':');
plot(thetane, fne, 'displayname', 'NE opt', 'color', colors(4,:), 'linestyle', ':');

% include start and end points
plot(thetape0(1), fpe0(1), 'o', 'color', colors(1,:), 'handlevisibility', 'off');
plot(thetape0(end), fpe0(end), 'x', 'color', colors(1,:), 'handlevisibility', 'off');
plot(thetane0(1), fne0(1), 'o', 'color', colors(2,:), 'handlevisibility', 'off');
plot(thetane0(end), fne0(end), 'x', 'color', colors(2,:), 'handlevisibility', 'off');
plot(thetape(1), fpe(1), 'o', 'color', colors(3,:), 'handlevisibility', 'off');
plot(thetape(end), fpe(end), 'x', 'color', colors(3,:), 'handlevisibility', 'off');
plot(thetane(1), fne(1), 'o', 'color', colors(4,:), 'handlevisibility', 'off');
plot(thetane(end), fne(end), 'x', 'color', colors(4,:), 'handlevisibility', 'off');

% include guest stoichiometries at start and end
plotGuestStoichiometries = true;
if plotGuestStoichiometries
    jsonInit = outputInit.jsonstruct;
    jsonOpt = outputOpt.jsonstruct;

    fsz = 14;
    style = {'linestyle', ':', 'handlevisibility', 'off'};
    left = {'horizontalalignment', 'left', 'fontsize', fsz, 'handlevisibility', 'off'};
    right = {'horizontalalignment', 'right', 'fontsize', fsz, 'handlevisibility', 'off'};

    theta0 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonInit.(pe).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [4, 5], 'color', colors(1,:), style{:});
    line([theta100, theta100], [4, 5], 'color', colors(1,:), style{:});
    text(theta0, 4.1, 'PE theta0 init', 'color', colors(1,:), left{:});
    text(theta100, 4.1, 'PE theta100 init', 'color', colors(1,:), left{:});

    theta0 = valsNotTruncated.(ne).guestStoichiometry0;
    %theta0 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonInit.(ne).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [0, 1], 'color', colors(2,:), style{:});
    line([theta100, theta100], [0, 1], 'color', colors(2,:), style{:});
    text(theta0, 0.1, 'NE theta0 init', 'color', colors(2,:), right{:});
    text(theta100, 0.1, 'NE theta100 init', 'color', colors(2,:), left{:});

    theta0 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonOpt.(pe).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [4, 5], 'color', colors(3,:), style{:});
    line([theta100, theta100], [4, 5], 'color', colors(3,:), style{:});
    text(theta0, 4.5, 'PE theta0 opt', 'color', colors(3,:), right{:});
    text(theta100, 4.5, 'PE theta100 opt', 'color', colors(3,:), right{:});

    theta0 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry0;
    theta100 = jsonOpt.(ne).(co).(am).(itf).guestStoichiometry100;
    line([theta0, theta0], [0, 1], 'color', colors(4,:), style{:});
    line([theta100, theta100], [0, 1], 'color', colors(4,:), style{:});
    text(theta0, 0.75, 'NE theta0 opt', 'color', colors(4,:), right{:});
    text(theta100, 0.75, 'NE theta100 opt', 'color', colors(4,:), left{:});

end

xlabel('theta')
title('cell balancing vs "lithiation"')
legend;

%%

return

%% Plot dQ/dV and Q(V)

t = expdata.time;
t = dataraw.time{1} * hour;
t = t(:);

% capacities: raw is a straight line (adjust sign)
capraw = -cumtrapz(dataraw.time{1}*hour, dataraw.current{1}); % straight, adjust sign

% processed: plot(cap) is not straight, but plot(expdata.time, cap) is straight
capexp = cumtrapz(expdata.time, ecs.expI*ones(size(expdata.time)));


% % extend time 10% before and after
% T = t(end) - t(1);
% t0 = t(1)-0.1*T;
% t1 = t(end)+0.1*T;
% t = linspace(t0, t1, numel(t));

[~, fpe0, fne0, thetape0, thetane0] = ecs.computeF(t, X0);
ocp0 = fpe0 - fne0;
[~, fpe, fne, thetape, thetane] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

colors = lines(4);

figure; hold on; grid on
plot(expdata.time/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
plot(t/hour, fpe0, 'displayname', 'pe 0', 'color', colors(1,:), 'linestyle', '--');
plot(t/hour, fne0, 'displayname', 'ne 0', 'color', colors(2,:), 'linestyle', '--');
plot(t/hour, ocp0, 'displayname', 'ocp 0', 'color', colors(3,:), 'linestyle', '--');
plot(t/hour, fpe, 'displayname', 'pe opt', 'color', colors(1,:));
plot(t/hour, fne, 'displayname', 'ne opt', 'color', colors(2,:));
plot(t/hour, ocp, 'displayname', 'ocp opt', 'color', colors(3,:));

xlabel 'Time  /  h';
ylabel 'Voltage  /  V';
legend('location', 'sw')

% plot derivatives
derivative = @(y,x) diff(y)./diff(x);
figure; hold on; grid on
% plot(t(2:end)/hour, derivative(fpe0, t), 'displayname', 'dV/dt pe 0', 'color', colors(1,:), 'linestyle', '--');
plot(t(2:end)/hour, derivative(fne0, t), 'displayname', 'dV/dt ne 0', 'color', colors(2,:), 'linestyle', '--');
%plot(t(2:end)/hour, derivative(ocp0, t), 'displayname', 'dV/dt ocp 0', 'color', colors(3,:), 'linestyle', '--');
%plot(t(2:end)/hour, derivative(fpe, t), 'displayname', 'dV/dt pe opt', 'color', colors(1,:));
plot(t(2:end)/hour, derivative(fne, t), 'displayname', 'dV/dt ne opt', 'color', colors(2,:));
%plot(t(2:end)/hour, derivative(ocp, t), 'displayname', 'dV/dt ocp opt', 'color', colors(3,:));
legend('location', 'nw')


% Plot vs capacity
t0 = t;
[~, fpe0, fne0] = ecs.computeF(t0, X0);
ocp0 = fpe0 - fne0;

[~, fpe, fne] = ecs.computeF(t, Xopt);
ocp = fpe - fne;

% figure; hold on; grid on
% plot(cap/hour, expdata.U, 'k--', 'displayname', 'Experiment 0.05 C');
% plot(cap/hour, fpe0, 'displayname', 'pe 0', 'color', colors(1,:), 'linestyle', '--');
% plot(cap/hour, fne0, 'displayname', 'ne 0', 'color', colors(2,:), 'linestyle', '--');
% plot(cap/hour, ocp0, 'displayname', 'ocp 0', 'color', colors(3,:), 'linestyle', '--');
% plot(cap/hour, fpe, 'displayname', 'pe opt', 'color', colors(1,:));
% plot(cap/hour, fne, 'displayname', 'ne opt', 'color', colors(2,:));
% plot(cap/hour, ocp, 'displayname', 'ocp opt', 'color', colors(3,:));

% xlabel 'Capacity  /  Ah';
% ylabel 'Voltage  /  V';
% legend('location', 'se')

% % plot dqdv
% figure; hold on; grid on
% %dqdvraw = gradient(capraw)./gradient(dataraw.voltage{1}, 2);
% dqdvexp = gradient(cap)./gradient(expdata.U);
% dqdv0 = gradient(cap)./gradient(ocp0);
% dqdv = gradient(cap)./gradient(ocp);
% %plot(dataraw.voltage{1}, dqdvraw, 'g', 'displayname', 'Experiment raw data');
% plot(expdata.U, -dqdvexp, 'k--', 'displayname', 'Experiment 0.05 C');
% plot(ocp0, -dqdv0, 'displayname', 'dQ/dV 0', 'color', colors(1,:), 'linestyle', '--');
% plot(ocp, -dqdv, 'displayname', 'dQ/dV opt', 'color', colors(1,:));
% minmax = @(x) [min(x), max(x)];
% disp(minmax(expdata.U));
% disp(minmax(ocp0));
% disp(minmax(ocp));
% disp(minmax(dqdvexp));
% disp(minmax(dqdv0));
% disp(minmax(dqdv));
% ylim([-2e2, 7e4])


% plot dqdv
qdvopts = {'nv', 500, ...
           'polyorder', 2, ...
           'difforder', 0, ...
           'x0', 0};
% [vgraw, qgraw, dqdvraw] = computedqdv(dataraw.time{1}*hour, dataraw.voltage{1}, abs(dataraw.current{1}), qdvopts{:});
% [vgexp, qgexp, dqdvexp] = computedqdv(expdata.time, expdata.U, ecs.expI*ones(size(expdata.time)), qdvopts{:});
% [vg0, qg0, dqdv0] = computedqdv(t0, ocp0, ecs.expI*ones(size(t0)), qdvopts{:});
% [vg, qg, dqdv] = computedqdv(t, ocp, ecs.expI*ones(size(t)), qdvopts{:});
% cap0 = cumtrapz(t0, ecs.expI*ones(size(t0)));
% cap = cumtrapz(t, ecs.expI*ones(size(t)));



% figure; hold on; grid on
% plot(vgexp)
% plot(vg0)
% plot(vg)

% figure; hold on; grid on
% plot(qgexp)
% plot(qg0)
% plot(qg)

% figure; hold on; grid on
% plot(vgexp, qgexp)
% plot(vg0, qg0)
% plot(vg, qg)

derivative = @(y,x) diff(y)./diff(x);

% plot q=q(v)
figure; hold on; grid on
plot(dataraw.voltage{1}, capraw, 'color', colors(3,:), 'displayname', 'Experiment raw data');
plot(vgraw, qgraw, 'o', 'color', colors(3,:), 'displayname', 'Experiment raw data (sgdf)');
plot(expdata.U, capexp, 'k--', 'displayname', 'Experiment');
plot(vgexp, qgexp, 'k--o', 'displayname', 'Experiment (sgdf)');
plot(ocp0, cap0, 'displayname', 'q(v) init', 'color', colors(1,:), 'linestyle', '--');
plot(vg0, qg0, 'o', 'displayname', 'q(v) init (sgdf)', 'color', colors(1,:), 'linestyle', '--');
plot(ocp, cap, 'displayname', 'q(v) opt', 'color', colors(2,:));
plot(vg, qg, 'o', 'displayname', 'q(v) opt (sgdf)', 'color', colors(2,:));
xlabel 'Voltage  /  V';
ylabel 'Capacity  /  As';
legend('location', 'sw')

figure; hold on; grid on
plot(vgraw, -gradient(qgraw)./gradient(vgraw), 'color', colors(3,:), 'displayname', 'Experiment raw data');
plot(vgexp, -gradient(qgexp)./gradient(vgexp), 'k--', 'displayname', 'Experiment 0.05 C');
plot(vg0, -gradient(qg0)./gradient(vg0), 'displayname', 'dQ/dV init', 'color', colors(1,:), 'linestyle', '--');
plot(vg, -gradient(qg)./gradient(vg), 'displayname', 'dQ/dV opt', 'color', colors(2,:));
title('dQ/dV computed from Q(V)');
legend('location', 'northwest')

figure; hold on; grid on
plot(vg, -gradient(qg)./gradient(vg), '.-', 'markersize', 15)
plot(vg(2:end), -derivative(qg, vg), 'o-')


% figure; hold on; grid on
% qmax = max(cap);
% plot(vgraw, dqdvraw, 'color', colors(3,:), 'displayname', 'Experiment raw data');
% plot(vgexp, dqdvexp, 'k--', 'displayname', 'Experiment 0.05 C');
% plot(vg0, dqdv0, 'displayname', 'dQ/dV 0', 'color', colors(1,:), 'linestyle', '--');
% plot(vg, dqdv, 'displayname', 'dQ/dV opt', 'color', colors(2,:));
% title(qdvopts);

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
