function HIT_Contour(cM2M2)
% Use 1.35 for the illustrative low-contact case.
if nargin<1, cM2M2=3.6; end
close all
set(0,'defaulttextinterpreter','latex','defaultaxesfontsize',16)
projectDir = setup_paths();
load(fullfile(projectDir,'results','DataFittingVariables.mat'))
final_param = final_param(~ismember({final_param.name}, {'rH','rM2'}));
idx = find(strcmp({final_param.name}, 'cM2'));
final_param(idx).name = 'cM2M2';
final_param(idx).value = 3.6;
n_params = numel(final_param);
p = struct();
for j = 1:n_params
    pname = final_param(j).name;
    p.(pname) = final_param(j).value;
end
p.qL = 0.8; p.qH = 0.8;
p.cM2M2 = cM2M2;
p.aM2 = 0.01;
vBase = struct();
vBase.L = (p.xiL +p.qL.*p.mu)./(p.xiL+p.mu);
vBase.H = (p.xiH +p.qH.*p.mu)./(p.xiH+p.mu);
[M1grid,M2grid] = meshgrid(0:0.01:1);
v = vBase;
v.M1 = M1grid;
v.M2 = M2grid;
rv = calculateR(p,v);
sFig{2}.name = 'Contour_R0_vm1vm2';
sFig{2}.figure = figure('Color', 'white', 'Position', [100 100 1400 650]);
vZero = struct();
vZero.L = 0;
vZero.H = 0;
vZero.M1 = M1grid;
vZero.M2 = M2grid;
rvZero = calculateR(p,vZero);
rvCompareMin = min([rvZero(:); rv(:)]);
rvCompareMax = max([rvZero(:); rv(:)]);
rvCompareLevels = linspace(0, ceil(rvCompareMax*4)/4, 25);
blackLevels = 0.5:0.3:rvCompareLevels(end);
blackLevels(blackLevels == 1) = [];
tCompare = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
comparisonData = {rvZero, rv};
comparisonTitles = {'(a) $v_l^* = 0$, $v_h^* = 0$', ...
    sprintf('(b) $v_l^* = %.2f$, $v_h^* = %.2f$', vBase.L, vBase.H)};
for panelIdx = 1:2
    ax = nexttile(tCompare);
    contourf(ax, M1grid, M2grid, comparisonData{panelIdx}, ...
        rvCompareLevels, 'LineStyle', 'none');
    hold(ax, 'on')
    [C,h] = contour(ax, M1grid, M2grid, comparisonData{panelIdx}, ...
        blackLevels, 'LineColor', [0.2 0.2 0.2], 'LineWidth', 2);
    clabel(C,h,'FontSize',12,'LabelSpacing',500)
    [D,g] = contour(ax, M1grid, M2grid, comparisonData{panelIdx}, ...
        [1 1], 'w-', 'LineWidth', 2.5);
    clabel(D,g,'FontSize',15,'Color','white','LabelSpacing',700)
    xlabel(ax, 'Vaccine Coverage of HetM ($v_{m1}^*$)')
    ylabel(ax, 'Vaccine Coverage of MSM ($v_{m2}^*$)')
    title(ax, comparisonTitles{panelIdx})
    clim(ax, [rvCompareLevels(1), rvCompareLevels(end)]);
    axis(ax, 'square')
    xlim(ax, [0 1]); ylim(ax, [0 1]);
    xticks(ax, 0:0.2:1); yticks(ax, 0:0.2:1);
    ax.Toolbar.Visible = 'off';
end
sgtitle(tCompare, '$\mathcal{R}_{v}$ over $(v_{m1}^*, v_{m2}^*)$ for two $(v_l^*, v_h^*)$ settings','interpreter','latex')
cbCompare = colorbar;
cbCompare.Layout.Tile = 'east';
cbCompare.Label.Interpreter = 'latex';
cbCompare.Label.String = '$\mathcal{R}_{v}$';
cbCompare.Ticks = 0:0.5:rvCompareLevels(end);
Hvals = [0 0.3 0.6 0.9];
M2vals = [0 0.3 0.6 0.9];
[Lgrid, M1gridLM] = meshgrid(0:0.01:1);
rvGrid = cell(numel(Hvals), numel(M2vals));
allRv = [];
for i = 1:numel(Hvals)
    for j = 1:numel(M2vals)
        vSweep = struct();
        vSweep.L = Lgrid;
        vSweep.H = Hvals(i);
        vSweep.M1 = M1gridLM;
        vSweep.M2 = M2vals(j);
        rvGrid{i,j} = calculateR(p,vSweep);
        allRv = [allRv; rvGrid{i,j}(:)];
    end
end
contourLevels = linspace(min(allRv), max(allRv), 12);
sFig{3}.name = 'Contour_R_vLvm1_by_vHvm2';
sFig{3}.figure = figure();
t = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',[100,100,1000,1000])
for i = 1:numel(Hvals)
    for j = 1:numel(M2vals)
        ax = nexttile;
        contourf(ax, Lgrid, M1gridLM, rvGrid{i,j}, contourLevels, 'LineStyle', 'none');
        clim(ax, [contourLevels(1), contourLevels(end)]);
        hold(ax, 'on')
        contour(ax, Lgrid, M1gridLM, rvGrid{i,j}, [1 1], 'w-', 'LineWidth', 1.5);
        xlim(ax, [0 1])
        ylim(ax, [0 1])
        xticks(ax, 0:0.3:1)
        yticks(ax, 0:0.3:1)
        axis(ax, 'square')
        if i == numel(Hvals)
            xlabel(ax, '$v_l^*$', 'FontSize', 20)
        end
        if j == 1
            ylabel(ax, '$v_{m1}^*$', 'FontSize', 20)
        end
    end
end
sgtitle(t, '$\mathcal{R}_{v}$ over $(v_l^*, v_{m1}^*)$ for fixed $v_h^*$ and $v_{m2}^*$','interpreter','latex')
cb = colorbar;
cb.Layout.Tile = 'east';
cb.FontSize = 18;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\mathcal{R}_{v}$';
cb.Label.FontSize = 22;
saveDir = fullfile(projectDir,'outputs','simulations',sprintf('contact_%g',cM2M2));
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
if cM2M2==3.6
    exportgraphics(sFig{2}.figure, [saveDir '/HIT_Contour_vm1vm2.png'], 'Resolution', 300);
end
exportgraphics(sFig{3}.figure, [saveDir '/HIT_Contour_vLvm1_grid.png'], 'Resolution', 300);
aM2Baseline = 0.01;
aM2PlausibleMax = 0.05;
msmExitE = p.sigM2 + p.mu;
msmExitI = p.psi_proM2 + p.psi_recM2 + p.mu;
msmExitP = p.g_proM2 + p.g_recM2 + p.mu;
msmInfectiousTime = p.etaM2/msmExitE ...
    + p.sigM2/(msmExitE*msmExitI) ...
    + p.thetaM2*p.sigM2*p.psi_proM2/(msmExitE*msmExitI*msmExitP);
Rmsm0 = p.betMM*p.cM2M2*msmInfectiousTime;
vM2Reference = (1 - 1/Rmsm0)/p.ev;
% The illustrative low-contact case has no positive within-MSM coverage bound.
if Rmsm0<=1, return; end
assert(vM2Reference >= 0 && vM2Reference <= 1, ...
    'The isolated MSM threshold is outside the plotted coverage range.');
fprintf('MSM vaccination coverage threshold: %.6f\n', vM2Reference);
aM2Values = unique([0:0.01:1, aM2Baseline]);
vM2Values = unique([0:0.01:1, vM2Reference]);
[aM2Grid, vM2Grid] = meshgrid(aM2Values, vM2Values);
pSweep = p;
pSweep.aM2 = aM2Grid;
vSweep = struct('L', 0.88, 'H', 0.88, 'M1', 0, 'M2', vM2Grid);
rvAM2M2 = calculateR(pSweep, vSweep);
aM2Limits = 0.5;
plotNamesAM2M2 = {'Contour_R_aM2_vM2', 'Contour_R_aM2_vM2_fullrange'};
fileNamesAM2M2 = {'HIT_Contour_aM2_vM2.png', ...
    'HIT_Contour_aM2_vM2_fullrange.png'};
for plotIdxAM2M2 = 1:numel(aM2Limits)
    figIdxAM2M2 = 3 + plotIdxAM2M2;
    aM2Max = aM2Limits(plotIdxAM2M2);
    inRangeAM2M2 = aM2Values <= aM2Max;
    rvDisplayedAM2M2 = rvAM2M2(:,inRangeAM2M2);
    rvLevelsAM2M2 = linspace(min(rvDisplayedAM2M2(:)), ...
        max(rvDisplayedAM2M2(:)), 12);
    sFig{figIdxAM2M2}.name = plotNamesAM2M2{plotIdxAM2M2};
    sFig{figIdxAM2M2}.figure = figure('Color', 'white', ...
        'Position', [100 100 800 800]);
    figAM2M2 = sFig{figIdxAM2M2}.figure;
    axAM2M2 = axes('Parent', figAM2M2);
    thresholdFontSize = 18;
    baselineFontSize = 16;
    plausibleFontSize = 14;
    if aM2Max == 0.5
        axAM2M2.FontSize = 20;
        axAM2M2.TitleFontSizeMultiplier = 1;
        thresholdFontSize = 22;
        baselineFontSize = 18;
        plausibleFontSize = 18;
    end
    contourf(axAM2M2, aM2Grid(:,inRangeAM2M2), ...
        vM2Grid(:,inRangeAM2M2), rvDisplayedAM2M2, ...
        rvLevelsAM2M2, 'LineStyle', 'none');
    clim(axAM2M2, [rvLevelsAM2M2(1), rvLevelsAM2M2(end)]);
    hold(axAM2M2, 'on')
    patch(axAM2M2, ...
        [0 aM2PlausibleMax aM2PlausibleMax 0], [0 0 1 1], ...
        'black', 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    xline(axAM2M2, aM2PlausibleMax, 'w:', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    C1 = contour(axAM2M2, aM2Grid(:,inRangeAM2M2), ...
        vM2Grid(:,inRangeAM2M2), rvAM2M2(:,inRangeAM2M2), ...
        [1 1], 'w-', 'LineWidth', 2, 'ShowText', 'off');
    xline(axAM2M2, aM2Baseline, 'y--', 'LineWidth', 2);
    yline(axAM2M2, vM2Reference, 'k--', 'LineWidth', 2);
    thresholdPoints = zeros(2,0);
    contourCol = 1;
    while contourCol < size(C1,2)
        nVertices = C1(2,contourCol);
        thresholdPoints = [thresholdPoints, ...
            C1(:,contourCol + (1:nVertices))];
        contourCol = contourCol + nVertices + 1;
    end
    labelX = 0.60*aM2Max;
    nearLabel = thresholdPoints(1,:) >= 0.46*aM2Max & ...
                thresholdPoints(1,:) <= 0.74*aM2Max;
    if any(nearLabel)
        labelY = max(thresholdPoints(2,nearLabel)) + 0.035;
        text(axAM2M2, labelX, labelY, '$\mathcal{R}_v=1$', ...
            'Interpreter', 'latex', 'Color', 'white', 'FontSize', thresholdFontSize, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    text(axAM2M2, 0.97*aM2Max, vM2Reference - 0.025, ...
        '$\mathcal{R}_{msm}=1$', 'Interpreter', 'latex', 'FontSize', thresholdFontSize, ...
        'Color', 'black', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top');
    xlabel(axAM2M2, '$a_{msm}$')
    ylabel(axAM2M2, '$v_{m2}^*$')
    title(axAM2M2, '$\mathcal{R}_v$ for $v_l^*=v_h^*=0.88$ and $v_{m_1}^*=0$')
    xlim(axAM2M2, [0 aM2Max])
    ylim(axAM2M2, [0 1])
    xticks(axAM2M2, linspace(0, aM2Max, 6))
    axis(axAM2M2, 'square')
    cbAM2M2 = colorbar(axAM2M2);
    cbAM2M2.Label.Interpreter = 'latex';
    cbAM2M2.Label.String = '$\mathcal{R}_v$';
    if aM2Max == 0.5
        cbAM2M2.FontSize = 20;
        cbAM2M2.Label.FontSize = 22;
    end
    cbAM2M2.Units = 'normalized';
    cbAM2M2.Position = [0.86 0.12 0.025 0.70];
    axAM2M2.Units = 'normalized';
    axAM2M2.Position = [0.13 0.12 0.70 0.70];
    if aM2Max == 0.5
        drawnow;
    end
    axPosAM2M2 = axAM2M2.Position;
    baselineX = axPosAM2M2(1) + axPosAM2M2(3)*aM2Baseline/aM2Max;
    axesTop = axPosAM2M2(2) + axPosAM2M2(4);
    annotation(figAM2M2, 'textarrow', ...
        [baselineX + 0.08, baselineX], ...
        axesTop - axPosAM2M2(4)*[0.06, 0.12], ...
        'String', 'Baseline: $a_{msm}=0.01$', 'Interpreter', 'latex', ...
        'FontSize', baselineFontSize, 'Color', 'yellow', 'LineWidth', 0.75, ...
        'HeadLength', 6, 'HeadWidth', 6, 'HorizontalAlignment', 'left');
    plausibleCenterX = axPosAM2M2(1) + ...
        axPosAM2M2(3)*(aM2PlausibleMax/2)/aM2Max;
    annotation(figAM2M2, 'textarrow', ...
        [axPosAM2M2(1) + 0.16*axPosAM2M2(3), plausibleCenterX], ...
        axPosAM2M2(2) + axPosAM2M2(4)*[0.46 0.40], ...
        'String', {'$a_{msm}\leq0.05$'}, ...
        'Interpreter', 'latex', 'FontSize', plausibleFontSize, ...
        'Color', 'white', 'LineWidth', 0.75, ...
        'HeadLength', 6, 'HeadWidth', 6, ...
        'HorizontalAlignment', 'left');
    exportgraphics(figAM2M2, ...
        fullfile(saveDir, fileNamesAM2M2{plotIdxAM2M2}), 'Resolution', 300);
end
end
