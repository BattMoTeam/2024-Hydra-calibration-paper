%% Export publication and supporting BattMo figures as .fig and .png files.

clear all
close all

mrstDebug(0);

set(0, 'defaultlinelinewidth', 2)
set(0, 'defaulttextfontsize', 14);
set(0, 'defaultaxesfontsize', 14);

am   = 'ActiveMaterial';
itf  = 'Interface';
pe   = 'PositiveElectrode';
ne   = 'NegativeElectrode';
co   = 'Coating';
sd   = 'SolidDiffusion';
ctrl = 'Control';

getTime = @(states) cellfun(@(s) s.time, states);
getE = @(states) cellfun(@(s) s.(ctrl).E, states);
getI = @(states) cellfun(@(s) s.(ctrl).I, states);

repoRoot = getHydra0Dir();
figuresDir = fullfile(repoRoot, 'figures');
supportingDir = fullfile(figuresDir, 'supporting');
ensureFolder(figuresDir);
ensureFolder(supportingDir);

saveddata = load(fullfile(repoRoot, 'raw-data', 'TE_1473.mat'));
dataraw = saveddata.experiment;

jsonstructEC = parseBattmoJson(fullfile(repoRoot, 'parameters', 'equilibrium-calibration-parameters.json'));
jsonstructHRC = parseBattmoJson(fullfile(repoRoot, 'parameters', 'high-rate-calibration-parameters.json'));
jsonstructHRC_1e14 = parseBattmoJson(fullfile(repoRoot, 'parameters', 'high-rate-calibration-parameters-1e-14-1e-14.json'));
jsonstructHRC_1e13 = parseBattmoJson(fullfile(repoRoot, 'parameters', 'high-rate-calibration-parameters-1e-13-1e-14.json'));

inputCap = struct('lowRateParams', jsonstructEC, ...
                  'include_current_collectors', true);
outputCap = runHydra(inputCap, 'runSimulation', false);
cap = computeCellCapacity(outputCap.model);

%% Figure 12. Cell balancing under equilibrium assumption

expdataLow = struct('time', dataraw.time{1} * hour, ...
                    'U'   , dataraw.voltage{1}    , ...
                    'I'   , abs(mean(dataraw.current{1})));

DRateLow = expdataLow.I / cap * hour;
outputLowInitial = runHydra(struct('DRate', DRateLow, ...
                                   'totalTime', expdataLow.time(end)), ...
                            'clearSimulation', true);
outputLowCalibrated = runHydra(struct('DRate', DRateLow, ...
                                      'totalTime', expdataLow.time(end), ...
                                      'lowRateParams', jsonstructEC), ...
                               'clearSimulation', true);
ecsInit = EquilibriumCalibrationSetup(outputLowInitial.model, expdataLow);
ecsOpt = EquilibriumCalibrationSetup(outputLowCalibrated.model, expdataLow);

t = expdataLow.time(:);
T = t(end) - t(1);
I = expdataLow.I;
tlong = linspace(t(1), t(end) + 0.4 * T, numel(t))';
q = cumtrapz(t, I * ones(size(t)));
qlong = cumtrapz(tlong, I * ones(size(tlong)));

[ocp0, fpe0, fne0] = ecsInit.computeF(tlong, ecsInit.X0);
[ocp, fpe, fne] = ecsOpt.computeF(t, ecsOpt.X0);

gne0 = outputLowInitial.jsonstruct.(ne).(co).(am).(itf).guestStoichiometry0;
gne100 = outputLowInitial.jsonstruct.(ne).(co).(am).(itf).guestStoichiometry100;
ocpminne = computeOCPanodeH0b(gne100);
ocpmaxne = computeOCPanodeH0b(gne0);
gpe0 = outputLowInitial.jsonstruct.(pe).(co).(am).(itf).guestStoichiometry0;
gpe100 = outputLowInitial.jsonstruct.(pe).(co).(am).(itf).guestStoichiometry100;
ocpminpe = computeOCPcathodeH0b(gpe0);

idx = find(fpe0 <= ocpminpe, 1, 'first');
qpe0cut = qlong(1:idx);
fpe0cut = fpe0(1:idx);
ocp0cut = ocp0(1:idx);

idx = find(fne0 >= ocpmaxne, 1, 'first');
qne0cut = qlong(1:idx);
fne0cut = fne0(1:idx);

tnelong = linspace(t(1) - 0.4 * T, tlong(end), numel(tlong))';
[~, ~, fnelong] = ecsOpt.computeF(tnelong, ecsOpt.X0);
qnelong = cumtrapz(tnelong, I * ones(size(tnelong))) - trapz([tnelong(1), t(1)], I * ones(2, 1));
idx = find(fnelong >= ocpmaxne, 1, 'first');
qnelong = qnelong(1:idx);
fnelong = fnelong(1:idx);
idx = find(fnelong >= 0, 1, 'first');
qnelong = qnelong(idx:end);
fnelong = fnelong(idx:end);

qscale = @(x) x / milli / outputLowInitial.jsonstruct.Geometry.faceArea * centi^2 / hour;
capExp_mAh_cm2 = qscale(q);
capPeInit_mAh_cm2 = qscale(qpe0cut);
capNeInit_mAh_cm2 = qscale(qne0cut);
capPeOpt_mAh_cm2 = qscale(q);
capNeOpt_mAh_cm2 = qscale(qnelong);

colors12 = lines(3);
fig12 = figure('Units', 'inches', 'Position', [0.2, 0.2, 7.2, 6.2]);
hold on
grid on
legend('Location', 'southwest')

plot(capNeInit_mAh_cm2, fne0cut, 'DisplayName', 'Graphite init', 'Color', colors12(1,:), 'LineStyle', '--');
plot(capPeInit_mAh_cm2, fpe0cut, 'DisplayName', 'LNMO init', 'Color', colors12(2,:), 'LineStyle', '--');
plot(capPeInit_mAh_cm2, ocp0cut, 'DisplayName', 'Full Cell init', 'Color', colors12(3,:), 'LineStyle', '--');
plot(capNeOpt_mAh_cm2, fnelong, 'DisplayName', 'Graphite opt', 'Color', colors12(1,:));
plot(capPeOpt_mAh_cm2, fpe, 'DisplayName', 'LNMO opt', 'Color', colors12(2,:));
plot(capPeOpt_mAh_cm2, ocp, 'DisplayName', 'Full Cell opt', 'Color', colors12(3,:));
plot(capExp_mAh_cm2, expdataLow.U, 'k:', 'DisplayName', 'Experiment 0.05 C');
xlabel('Capacity / mAh cm^{-2}')
ylabel('Voltage / V')
title('Figure 12. Cell balancing under equilibrium assumption')
axis tight
breakyaxis([1.5, 3]);

saveFigureSet(fig12, fullfile(figuresDir, 'figure-12-cell-balancing-under-equilibrium-assumption'));
saveJson( ...
    fullfile(figuresDir, 'figure-12-cell-balancing-under-equilibrium-assumption.json'), ...
    struct( ...
    'title', 'Figure 12. Cell balancing under equilibrium assumption', ...
    'capacity_unit', 'mAh cm^-2', ...
    'experiment', struct( ...
    'capacity_mAh_cm2', capExp_mAh_cm2(:)', ...
    'voltage_v', expdataLow.U(:)'), ...
    'graphite_init', struct( ...
    'capacity_mAh_cm2', capNeInit_mAh_cm2(:)', ...
    'voltage_v', fne0cut(:)'), ...
    'lnmo_init', struct( ...
    'capacity_mAh_cm2', capPeInit_mAh_cm2(:)', ...
    'voltage_v', fpe0cut(:)'), ...
    'full_cell_init', struct( ...
    'capacity_mAh_cm2', capPeInit_mAh_cm2(:)', ...
    'voltage_v', ocp0cut(:)'), ...
    'graphite_opt', struct( ...
    'capacity_mAh_cm2', capNeOpt_mAh_cm2(:)', ...
    'voltage_v', fnelong(:)'), ...
    'lnmo_opt', struct( ...
    'capacity_mAh_cm2', capPeOpt_mAh_cm2(:)', ...
    'voltage_v', fpe(:)'), ...
    'full_cell_opt', struct( ...
    'capacity_mAh_cm2', capPeOpt_mAh_cm2(:)', ...
    'voltage_v', ocp(:)')) ...
        );

%% Figure 13. Results after high-rate calibration at 2C using two different initial values...

khigh = numel(dataraw.time);
expdataHigh = struct('time', dataraw.time{khigh} * hour, ...
                     'U'   , dataraw.voltage{khigh}    , ...
                     'I'   , abs(mean(dataraw.current{khigh})));

DRateHigh = expdataHigh.I / cap * hour;
commonHighInput = struct('DRate'                         , DRateHigh           , ...
                         'totalTime'                     , expdataHigh.time(end), ...
                         'lowRateParams'                 , jsonstructEC        , ...
                         'useRegionBruggemanCoefficients', true                , ...
                         'include_current_collectors'    , true);

inputHighInitial1 = commonHighInput;
inputHighInitial1.Dne = 1e-14;
inputHighInitial1.Dpe = 1e-14;
inputHighInitial2 = commonHighInput;
inputHighInitial2.Dne = 1e-13;
inputHighInitial2.Dpe = 1e-14;
inputHighCal1 = commonHighInput;
inputHighCal1.highRateParams = jsonstructHRC_1e14;
inputHighCal2 = commonHighInput;
inputHighCal2.highRateParams = jsonstructHRC_1e13;

outputHighInitial1 = runHydra(inputHighInitial1, 'clearSimulation', true);
outputHighInitial2 = runHydra(inputHighInitial2, 'clearSimulation', true);
outputHighCal1 = runHydra(inputHighCal1, 'clearSimulation', true);
outputHighCal2 = runHydra(inputHighCal2, 'clearSimulation', true);

fig13 = figure('Units', 'inches', 'Position', [0.2, 0.2, 8.2, 5.8]);
hold on
grid on

plot(expdataHigh.time / hour, expdataHigh.U, 'k--', 'DisplayName', 'Experiment 2C');
plot(getTime(outputHighInitial1.states) / hour, getE(outputHighInitial1.states), 'Color', [0.75, 0.45, 0.1], 'DisplayName', 'Initial D_{NE}=10^{-14}');
plot(getTime(outputHighInitial2.states) / hour, getE(outputHighInitial2.states), 'Color', [0.5, 0.25, 0.7], 'DisplayName', 'Initial D_{NE}=10^{-13}');
plot(getTime(outputHighCal1.states) / hour, getE(outputHighCal1.states), '-', 'Color', [0.1, 0.45, 0.85], 'DisplayName', 'Calibrated from D_{NE}=10^{-14}');
plot(getTime(outputHighCal2.states) / hour, getE(outputHighCal2.states), '--', 'Color', [0.85, 0.2, 0.2], 'DisplayName', 'Calibrated from D_{NE}=10^{-13}');
xlabel('Time / h')
ylabel('Voltage / V')
title('Figure 13. Results after high-rate calibration at 2C')
legend('Location', 'southwest')
ylim([3.2, 4.95])
saveFigureSet(fig13, fullfile(figuresDir, 'figure-13-high-rate-calibration-at-2C'));
saveJson( ...
    fullfile(figuresDir, 'figure-13-high-rate-calibration-at-2C.json'), ...
    struct( ...
    'title', 'Figure 13. Results after high-rate calibration at 2C using two different initial values for the diffusion coefficient of the negative electrode. The calibrated curves lie on top of each other', ...
    'experiment', struct('time_h', expdataHigh.time(:)' / hour, 'voltage_v', expdataHigh.U(:)'), ...
    'initial_guess_dne_1e_14', struct('time_h', getTime(outputHighInitial1.states) / hour, 'voltage_v', getE(outputHighInitial1.states)), ...
    'initial_guess_dne_1e_13', struct('time_h', getTime(outputHighInitial2.states) / hour, 'voltage_v', getE(outputHighInitial2.states)), ...
    'calibrated_from_dne_1e_14', struct('time_h', getTime(outputHighCal1.states) / hour, 'voltage_v', getE(outputHighCal1.states)), ...
    'calibrated_from_dne_1e_13', struct('time_h', getTime(outputHighCal2.states) / hour, 'voltage_v', getE(outputHighCal2.states))) ...
        );

%% Figure 14. Experimental voltages and P2D model results over different discharge rates

cases = cell(numel(dataraw.time), 1);
colors = lines(numel(dataraw.time));

for k = 1:numel(dataraw.time)
    expdata = struct('time', dataraw.time{k} * hour, ...
                     'U'   , dataraw.voltage{k}    , ...
                     'I'   , abs(mean(dataraw.current{k})));

    DRate = expdata.I / cap * hour;
    output = runHydra(struct('DRate'                         , DRate             , ...
                             'totalTime'                     , expdata.time(end) , ...
                             'lowRateParams'                 , jsonstructEC      , ...
                             'highRateParams'                , jsonstructHRC     , ...
                             'useRegionBruggemanCoefficients', true              , ...
                             'include_current_collectors'    , true), ...
                      'clearSimulation', true);

    cases{k} = buildCaseStruct(sprintf('Discharge rate %d', k), expdata, output, DRate, colors(k,:), ne, pe, co, am, sd, itf, ctrl);
end

fig14 = figure('Units', 'inches', 'Position', [0.2, 0.2, 8.2, 5.8]);
hold on
grid on

hp2d = gobjects(numel(cases), 1);
for k = 1:numel(cases)
    caseData = cases{k};
    plot(caseData.exp_capacity_ah, caseData.experimental_voltage_v, '--', 'Color', caseData.color);
    hp2d(k) = plot(caseData.sim_capacity_ah, caseData.sim_voltage_v, '-', 'Color', caseData.color);
end

xlabel('Capacity / Ah')
ylabel('Voltage / V')
title('Figure 14. Experimental voltages and P2D model results over different discharge rates')
ylim([3.45, 4.9])
axis tight

hp(1) = plot(nan, nan, 'k--'); %#ok<AGROW>
hp(2) = plot(nan, nan, 'k-'); %#ok<AGROW>
legend(gca, hp, {'Experiment', 'P2D'}, 'Location', 'northwest')

legtxt = arrayfun(@(caseData) sprintf('%1.2gC', caseData.drate), [cases{:}], 'UniformOutput', false);
ax = axes('Position', get(gca, 'Position'), 'Visible', 'off');
legend(ax, hp2d, legtxt, 'Location', 'southwest')

saveFigureSet(fig14, fullfile(figuresDir, 'figure-14-experimental-voltages-and-p2d-results'));
saveJson( ...
    fullfile(figuresDir, 'figure-14-experimental-voltages-and-p2d-results.json'), ...
    struct( ...
    'title', 'Figure 14. Experimental voltages and P2D model results over different discharge rates.', ...
    'cases', {cellfun(@caseToJsonStruct, cases, 'UniformOutput', false)}) ...
        );

%% Supporting voltage curves and state dashboards for each validation run

for k = 1:numel(cases)
    caseData = cases{k};
    caseSlug = slugify(caseData.case_name);
    caseDir = fullfile(supportingDir, caseSlug);
    ensureFolder(caseDir);

    figVoltage = figure('Units', 'inches', 'Position', [0.2, 0.2, 8.2, 5.4]);
    hold on
    grid on
    plot(caseData.exp_capacity_ah, caseData.experimental_voltage_v, 'k--', 'DisplayName', 'Experiment');
    plot(caseData.sim_capacity_ah, caseData.sim_voltage_v, 'Color', caseData.color, 'DisplayName', 'BattMo');
    xlabel('Capacity / Ah')
    ylabel('Voltage / V')
    title(sprintf('%s voltage curve (%1.3g A, %1.2gC)', caseData.case_name, caseData.current_a, caseData.drate))
    legend('Location', 'southwest')
    ylim([3.35, 4.95])
    saveFigureSet(figVoltage, fullfile(caseDir, sprintf('%s-voltage', caseSlug)));

    figDashboard = figure('Units', 'inches', 'Position', [0.2, 0.2, 14, 8.5]);
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    plotStateTile(nexttile, caseData.time_h, caseData.x_elyte_um, caseData.elyte_c, 'Electrolyte concentration / mol m^{-3}');
    plotStateTile(nexttile, caseData.time_h, caseData.x_elyte_um, caseData.elyte_phi, 'Electrolyte potential / V');
    plotStateTile(nexttile, caseData.time_h, caseData.x_ne_um, caseData.ne_phi, 'Negative electrode potential / V');
    plotStateTile(nexttile, caseData.time_h, caseData.x_pe_um, caseData.pe_phi, 'Positive electrode potential / V');
    plotStateTile(nexttile, caseData.time_h, caseData.x_ne_um, caseData.ne_theta, 'Negative particle surface stoichiometry / -');
    plotStateTile(nexttile, caseData.time_h, caseData.x_pe_um, caseData.pe_theta, 'Positive particle surface stoichiometry / -');
    sgtitle(sprintf('%s state variable contour dashboard', caseData.case_name))
    saveFigureSet(figDashboard, fullfile(caseDir, sprintf('%s-state-dashboard', caseSlug)));
end

saveJson( ...
    fullfile(supportingDir, 'battmo-validation-states.json'), ...
    struct( ...
    'title', 'BattMo validation states for supporting interactive documentation', ...
    'cases', {cellfun(@caseToJsonStruct, cases, 'UniformOutput', false)}) ...
        );

fprintf('Wrote publication figures and supporting dashboards to %s\n', figuresDir);


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
