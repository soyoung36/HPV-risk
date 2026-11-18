function AvertedCancers()
% Compare group-targeted vaccination and spillover over 2021-2051.
close all
set(0,'defaulttextinterpreter','latex','defaultaxesfontsize',16)
projectDir = setup_paths();
load(fullfile(projectDir,'results','DataFittingVariables.mat'))
final_param = final_param(~ismember({final_param.name}, {'rM2','rH'}));
idx = find(strcmp({final_param.name}, 'cM2'));
final_param(idx).name = 'cM2M2';
final_param(idx).value = 3.6;
historicalStartYear = data.time.all(1);
switchTime = 2007 - historicalStartYear;
projectionStartTime = 2021 - historicalStartYear;
p = struct();
for i = 1:numel(final_param)
    p.(final_param(i).name) = final_param(i).value;
end
historicalInitial = [initial.value]';
historicalInitial([1 9]) = historicalInitial([1 9]) + historicalInitial([2 10]);
historicalInitial([2 10]) = 0;
pPre = p;
pPre.qL = 0; pPre.qH = 0;
pPre.xiL = 0; pPre.xiH = 0;
[tPre, yPre] = ode45(@(t,y) ODE_Historical(t,y,pPre), ...
    0:switchTime, historicalInitial);
pPost = p;
pPost.qL = 0.05; pPost.qH = 0.05;
pPost.xiL = 0.01; pPost.xiH = 0.01;
[tPost, yPost] = ode45(@(t,y) ODE_Historical(t,y,pPost), ...
    switchTime:projectionStartTime, yPre(end,:)');
tHistorical = [tPre; tPost(2:end)];
ySim = [yPre; yPost(2:end,:)];
temp = ySim(end,:);
initialvalue_2021 = [temp(1:end-3), temp(end-2)*0.8, temp(end-2)*0.2, temp(end-1)*0.8, temp(end-1)*0.2,...
    temp(end)*0.5, temp(end)*0.5]';
p.qL = 0; p.qH = 0; p.qM1 = 0; p.qM2 = 0;
p.xiL = 0; p.xiH = 0; p.xiM1 = 0; p.xiM2 = 0;
tspan = 0:30;
[~, ySim] = ode45(@(t,y) ODE_4groupIncidence(t,y,p), tspan, initialvalue_2021);
sim_fields = {'CCl','CCh','CDl','CDh','CCm1','CCm2','CDf'};
for f = 1:6
    ctrl.(sim_fields{f}) = ySim(:,end-6+f);
end
ctrl.CDf = ctrl.CDl + ctrl.CDh;
p0 = p;
vac_fields = {'qL','qH','qM1','qM2'};
for f = 1:7
    averted.(sim_fields{f}) = [];
    averted_ratio.(sim_fields{f}) = [];
end
for j = 1:4
    p = p0;
    p.(vac_fields{j}) = 0.8;
    [~, ySim] = ode45(@(t,y) ODE_4groupIncidence(t,y,p), tspan, initialvalue_2021);
    for f = 1:6
        averted.(sim_fields{f}) = [averted.(sim_fields{f}), ctrl.(sim_fields{f}) - ySim(:,end-6+f)];
        averted_ratio.(sim_fields{f}) = [averted_ratio.(sim_fields{f}), averted.(sim_fields{f})(:,end)./ctrl.(sim_fields{f})];
    end
    averted.CDf = [averted.CDf, ctrl.CDl+ctrl.CDh - ySim(:,end-3)-ySim(:,end-2)];
    averted_ratio.CDf = [averted_ratio.CDf, averted.CDf(:,end)./ctrl.CDf];
end
aFig{1}.name = 'AvertedCancerCases';
aFig{1}.figure = figure();
set(gcf,'Position',[100,100,1100,700])
titles = {'(a) Averted Cancer Cases in Low-risk Females','(b) Averted Cancer Cases in High-risk Females','(c) Averted Cancer Cases in Heterosexual Males','(d) Averted Cancer Cases in MSM'};
fig_fields = {'CCl','CCh','CCm1','CCm2'};
cancerAxes = gobjects(4,1);
for j = 1:4
    cancerAxes(j) = subplot(2,2,j);
    f1 = plot(tspan+2021, averted.(fig_fields{j})', 'LineWidth', 3.5);
    title(titles{j});
    xlabel('Time (years)');
    ylabel('Averted Cancer Cases');
    set(gca, 'YDir','reverse')
    grid on;
    if j == 1
        cancerLegendLines = f1;
    end
end
cancerLegend = legend(cancerAxes(1), cancerLegendLines, ...
    {'Vaccination of low-risk females', 'Vaccination of high-risk females', ...
     'Vaccination of heterosexual males', 'Vaccination of MSM'}, ...
    'NumColumns', 2, 'FontSize', 13, 'Interpreter', 'none', ...
    'Location', 'none', 'Box', 'off');
cancerLegend.Units = 'pixels';
drawnow;
figurePosition = aFig{1}.figure.Position;
legendPosition = cancerLegend.Position;
legendPosition(1) = (figurePosition(3) - legendPosition(3))/2;
legendPosition(2) = 6;
cancerLegend.Position = legendPosition;
set(cancerAxes, 'Units', 'pixels');
cancerPanelPositions = vertcat(cancerAxes.Position);
extraHeight = 0;
for j = 3:4
    labelInsets = cancerAxes(j).TightInset;
    requiredBottom = legendPosition(2) + legendPosition(4) + labelInsets(2) + 8;
    extraHeight = max(extraHeight, requiredBottom - cancerPanelPositions(j,2));
end
extraHeight = ceil(extraHeight);
figurePosition(2) = max(0, figurePosition(2) - extraHeight);
figurePosition(4) = figurePosition(4) + extraHeight;
aFig{1}.figure.Position = figurePosition;
for j = 1:4
    panelPosition = cancerPanelPositions(j,:);
    panelPosition(2) = panelPosition(2) + extraHeight;
    cancerAxes(j).Position = panelPosition;
end
aFig{2}.name = 'AvertedCancerDeaths';
aFig{2}.figure = figure();
set(gcf,'Position',[0,50,1700,400])
titles_deaths = {'(a) Averted Cancer Deaths in Low-risk Females','(b) Averted Cancer Deaths in High-risk Females','(c) Averted Cancer Deaths in Females'};
fig_fields_death = {'CDl','CDh','CDf'};
deathAxes = gobjects(3,1);
for k = 1:3
    deathAxes(k) = subplot(1,3,k);
    f2 = plot(tspan+2021, averted.(fig_fields_death{k})', 'LineWidth', 3.5);
    title(titles_deaths{k});
    xlabel('Time (years)');
    ylabel('Averted Cancer Deaths');
    set(gca, 'YDir','reverse')
    grid on;
    if k == 1
        deathLegendLines = f2;
    end
end
deathLegend = legend(deathAxes(1), deathLegendLines, ...
    {'Vaccination of low-risk females', 'Vaccination of high-risk females', ...
     'Vaccination of heterosexual males', 'Vaccination of MSM'}, ...
    'NumColumns', 2, 'FontSize', 13, 'Interpreter', 'none', ...
    'Location', 'none', 'Box', 'off');
deathLegend.Units = 'pixels';
drawnow;
deathFigurePosition = aFig{2}.figure.Position;
deathLegendPosition = deathLegend.Position;
deathLegendPosition(1) = (deathFigurePosition(3) - deathLegendPosition(3))/2;
deathLegendPosition(2) = 6;
deathLegend.Position = deathLegendPosition;
set(deathAxes, 'Units', 'pixels');
deathPanelPositions = vertcat(deathAxes.Position);
deathExtraHeight = 0;
for k = 1:3
    deathLabelInsets = deathAxes(k).TightInset;
    deathRequiredBottom = deathLegendPosition(2) + deathLegendPosition(4) ...
        + deathLabelInsets(2) + 8;
    deathExtraHeight = max(deathExtraHeight, ...
        deathRequiredBottom - deathPanelPositions(k,2));
end
deathExtraHeight = ceil(deathExtraHeight);
deathFigurePosition(2) = max(0, deathFigurePosition(2) - deathExtraHeight);
deathFigurePosition(4) = deathFigurePosition(4) + deathExtraHeight;
aFig{2}.figure.Position = deathFigurePosition;
for k = 1:3
    deathPanelPosition = deathPanelPositions(k,:);
    deathPanelPosition(2) = deathPanelPosition(2) + deathExtraHeight;
    deathAxes(k).Position = deathPanelPosition;
end
set(deathAxes, 'Units', 'normalized');
deathLegend.Units = 'normalized';
aM2Values = [0.005, 0.01, 0.02, 0.05];
projectionTime = [0 30];
sweepOptions = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
initialSweep = initialvalue_2021;
initialSweep(33:38) = 0;
cancerIdx = [33 34 37 38];
groupLabels = {'Low-risk females', 'High-risk females', ...
    'Heterosexual males', 'MSM'};
strategyLabels = {'Vaccination of low-risk females', ...
    'Vaccination of high-risk females', ...
    'Vaccination of heterosexual males', 'Vaccination of MSM'};
percentAverted = nan(numel(aM2Values), 4, 4);
controlCases = zeros(numel(aM2Values), 4);
vaccinatedCases = zeros(numel(aM2Values), 4, 4);
for aIdx = 1:numel(aM2Values)
    pControl = p0;
    pControl.aM2 = aM2Values(aIdx);
    [~, yControl] = ode45(@(t,y) ODE_4groupIncidence(t,y,pControl), ...
        projectionTime, initialSweep, sweepOptions);
    controlCases(aIdx,:) = yControl(end,cancerIdx);
    for strategyIdx = 1:4
        pVaccinated = pControl;
        pVaccinated.(vac_fields{strategyIdx}) = 0.8;
        [~, yVaccinated] = ode45( ...
            @(t,y) ODE_4groupIncidence(t,y,pVaccinated), ...
            projectionTime, initialSweep, sweepOptions);
        vaccinatedCases(aIdx,:,strategyIdx) = yVaccinated(end,cancerIdx);
        validGroups = controlCases(aIdx,:) > 0;
        percentAverted(aIdx,validGroups,strategyIdx) = 100 * ...
            (controlCases(aIdx,validGroups) ...
            - vaccinatedCases(aIdx,validGroups,strategyIdx)) ...
            ./ controlCases(aIdx,validGroups);
    end
end
aFig{3}.name = 'AvertedCancerCases_aM2_Comparison';
aFig{3}.figure = figure('Color', 'white', ...
    'Position', [100 100 1700 650], ...
    'DefaultAxesFontSize', 24, ...
    'DefaultAxesLabelFontSizeMultiplier', 1, ...
    'DefaultAxesTitleFontSizeMultiplier', 1, ...
    'DefaultLegendFontSize', 22);
tAverted = tiledlayout(aFig{3}.figure, 1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(tAverted, 1);
plot(aM2Values, percentAverted(:,1:3,4), '-o', ...
    'LineWidth', 3.5, 'MarkerSize', 7);
xlabel('$a_{msm}$')
ylabel('Cumulative cancer cases averted (\%)')
title('(a) MSM-only vaccination during 2021-2051')
xticks(aM2Values)
xticklabels({'0.005','0.01','0.02','0.05'})
xlim([aM2Values(1), aM2Values(end)])
legend(groupLabels(1:3), 'Location', 'northwest', 'Interpreter', 'none');
grid on
box on
nexttile(tAverted, 2);
plot(aM2Values, squeeze(percentAverted(:,4,:)), '-o', ...
    'LineWidth', 3.5, 'MarkerSize', 7);
xlabel('$a_{msm}$')
ylabel('Cumulative MSM cancer cases averted (\%)')
title('(b) Single group vaccination during 2021-2051')
xticks(aM2Values)
xticklabels({'0.005','0.01','0.02','0.05'})
xlim([aM2Values(1), aM2Values(end)])
legend(strategyLabels, 'Location', 'best', 'Interpreter', 'none');
grid on
box on
saveDir = fullfile(projectDir,'outputs','simulations');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
for afIdx=1:numel(aFig)
    exportgraphics(aFig{afIdx}.figure, [saveDir '/' aFig{afIdx}.name '.png'], 'Resolution', 300);
end
end
