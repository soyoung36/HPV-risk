function SensitivityAnalysis(resp_ftn,fresh,seed)
% fresh=true repeats screening; nonbaseline seeds with fresh=false reuse saved selections.
if nargin<1, resp_ftn='R'; end
if nargin<2, fresh=false; end
if nargin<3, seed=1; end
validateattributes(seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
assert(ismember(resp_ftn,{'R','Cf','Cm','Cmsm','all'}),'Unknown response.');
close all;
set(0, 'defaultaxesfontsize',25,'defaulttextinterpreter','latex')
projectDir = setup_paths();
if ~fresh && (seed~=1 || strcmp(resp_ftn,'all'))
    repeatPRCC(projectDir,resp_ftn,seed);
    return
end
assert(~strcmp(resp_ftn,'all'),'Use fresh=false with all to reuse the saved screening.');
load(fullfile(projectDir,'results','DataFittingVariables.mat'))
sen_param = final_param(1:end-2);
for i = 5:8
    sen_param(i).value = 0.5;
    sen_param(i+25).value = 0.01;
end
sen_param = sen_param(~ismember({sen_param.name}, {'rH','rM2'}));
idx = find(strcmp({sen_param.name}, 'cM2'));
sen_param(idx).name = 'cM2M2';
sen_param(idx).value = 3.6;
percent = 0.2;
for j = 1:numel(sen_param)
    sen_param(j).range = [sen_param(j).value * (1 - percent), sen_param(j).value * (1 + percent)];
end
sen_param(strcmp({sen_param.name},'ev')).range(2) = 1;
plabel = {'$\pi_l$','$\pi_h$','$\pi_{m1}$','$\pi_{m2}$','$q_l$','$q_h$','$q_{m1}$','$q_{m2}$',...
    '$\zeta_l$','$\zeta_h$','$\zeta_{m1}$','$\zeta_{m2}$','$\beta_{mf}$','$\beta_{fm}$','$\beta_{mm}$',...
    '$c_{hetf_l}$','$c_{msm}^{msm}$','$a_{msm}$','$a_{hetf_h}$','$\eta_l$','$\eta_h$','$\eta_{m1}$','$\eta_{m2}$',...
    '$\theta_l$','$\theta_h$','$\theta_{m1}$','$\theta_{m2}$','$\xi_l$','$\xi_h$','$\xi_{m1}$','$\xi_{m2}$',...
    '$\mu$','$\varepsilon_v$','$\sigma_l$','$\sigma_h$','$\sigma_{m1}$','$\sigma_{m2}$',...
    '$(1-b_l)\psi_l$','$(1-b_h)\psi_h$','$(1-b_{m1})\psi_{m1}$','$(1-b_{m2})\psi_{m2}$',...
    '$b_l\psi_l$','$b_h\psi_h$','$b_{m1}\psi_{m1}$','$b_{m2}\psi_{m2}$',...
    '$(1-d_l)g_l$','$(1-d_h)g_h$','$(1-d_{m1})g_{m1}$','$(1-d_{m2})g_{m2}$','$d_lg_l$','$d_hg_h$','$d_{m1}g_{m1}$','$d_{m2}g_{m2}$',...
    '$w_l$','$w_h$','$w_{m1}$','$w_{m2}$','$\delta_l$','$\delta_h$','$\delta_{m1}$','$\delta_{m2}$'};
response = struct();
t_end = 50;
switch resp_ftn
    case 'R'
        response.label = '$R_v$';
    case 'Cf'
        response.label = '$C_f$';
    case 'Cm'
        response.label = '$C_m$';
    case 'Cmsm'
        response.label = '$C_{msm}$';
end
if fresh
rng(seed)
nLHS = 1000;
nPar = numel(sen_param);
for i = 1:nPar
    sen_param(i).LHS = LHS_Call(sen_param(i).range(1), sen_param(i).value, sen_param(i).range(2),0,nLHS,'unif');
end
resp_vals = cell(nPar, 1);
idx_in = [];
idx_influ = [];
influence_threshold = 0.05;
p_sen = struct();
for j = 1:nPar
    pname = sen_param(j).name;
    p_sen.(pname) = sen_param(j).value;
end
disp('Parameters')
switch resp_ftn
    case 'R'
    for i = 1:nPar
        p = p_sen;
        pname = sen_param(i).name;
        p.(pname) = sen_param(i).LHS;
        resp_vals{i} = calculateR(p);
        resp = resp_vals{i};
        response_is_constant = isscalar(resp);
        if response_is_constant
            fprintf('%s,',sen_param(i).name)
        else idx_in(end+1) = i;
        end
        if var(resp) >= 1e-6
            idx_influ(end+1) = i;
        end
    end
    disp('are not included in R')
    disp('Parameters')
    for i = 1:length(idx_influ)
        j = idx_influ(i);
        fprintf('%s, ',sen_param(j).name)
    end
    disp('are influential on R')
    case {'Cf','Cm','Cmsm'}
    for i = 1:nPar
        Cvals = zeros(nLHS,1);
        pname = sen_param(i).name;
        pLHS = sen_param(i).LHS;
        p = p_sen;
        for k = 1:nLHS
            p.(pname) = pLHS(k);
            C_end = calculateC(initial,p,t_end);
            switch resp_ftn
                case 'Cf'
                    Cvals(k) = C_end.Cf;
                case 'Cm'
                    Cvals(k) = C_end.Cm;
                case 'Cmsm'
                    Cvals(k) = C_end.Cmsm;
            end
        end
        resp_vals{i} = Cvals;
        if (max(Cvals)-min(Cvals))/mean(Cvals) >= influence_threshold
            idx_influ(end+1) = i;
            fprintf('%s, ',sen_param(i).name)
        end
    end
    disp(['are influential on the cancer at t = ',num2str(t_end)])
end
else
    load(fullfile(projectDir,'results',['Sensitivity_' resp_ftn '.mat']));
end
close all
fprintf('Generating scatter plots...\n');
nplot = 16;
nrow = 2;
ncol = 8;
scatter_option = 'influ';
switch scatter_option
    case 'all'
        nScatter = nPar;
    case 'influ'
        nScatter = numel(idx_influ);
end
nfig = ceil(nScatter / nplot);
aFig = cell(nfig, 1);
for i = 1:nScatter
    figIdx = ceil(i / nplot);
    subplotIdx = mod(i - 1, nplot) + 1;
    if subplotIdx == 1
        aFig{figIdx}.name = ['ScatterPlots_',resp_ftn,num2str(figIdx)];
        aFig{figIdx}.figure = figure(figIdx);
        set(gcf,'Position',[100,100,1600,500])
        if figIdx == 1
            if resp_ftn == 'R'
                sgtitle(['One-at-a-time Response Curves of ',response.label], 'FontSize', 30);
            else
                sgtitle(['One-at-a-time Response Curves of ',strrep(response.label,'C_','J_'),' at t = ',num2str(t_end)],'FontSize',30);
            end
        end
    end
    switch scatter_option
        case 'all'
            pIdx = i;
        case 'influ'
            pIdx = idx_influ(i);
    end
    subplot(nrow, ncol, subplotIdx);
    scatter(sen_param(pIdx).LHS, resp_vals{pIdx}, 'filled');
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    ax.XAxis.Exponent = 0;
    xlabel([plabel(pIdx)]);
    grid on;
    if strcmp(resp_ftn, 'R')
        crowdedParameters = {'betMM','thetaM2','xiM2','mu','psi_proM2','g_recM2'};
    else
        crowdedParameters = {'piL','piH','piM1','piM2','betFM','betMM', ...
            'aM2','thetaL','thetaH','thetaM1','xiL','xiM1','xiM2','mu', ...
            'g_recL','g_recH','psi_proM1','psi_proM2'};
    end
    if ismember(sen_param(pIdx).name, crowdedParameters)
        ticks = ax.XTick;
        ax.XTick = ticks([1,end]);
        ax.XTickLabelRotation = 0;
        if strcmp(sen_param(pIdx).name, 'piM2')
            ax.XTick = 1000 * [ceil(ax.XLim(1)/1000), floor(ax.XLim(2)/1000)];
        end
        if startsWith(sen_param(pIdx).name, 'pi')
            ax.XTickLabel = arrayfun(@(x) regexprep(sprintf('%.2g',x), ...
                'e\+0*','e'), ax.XTick, 'UniformOutput', false);
        end
    end
end
fprintf('Finished generating %d scatter plots.\n', nPar);
saveDir = fullfile(projectDir,'outputs','sensitivity',['scatter_' resp_ftn]);
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
for f = 1:nfig
    exportgraphics(aFig{f}.figure, [saveDir '/' aFig{f}.name '.png'], 'Resolution', 300);
    exportgraphics(aFig{f}.figure, ...
        fullfile(saveDir, ['OAT_Response_Curves_' resp_ftn '.pdf']), ...
        'ContentType', 'vector', 'Append', f > 1);
end
if fresh
LHSmat = zeros(nLHS,nPar);
for j = 1:nPar
    LHSmat(:,j) = sen_param(j).LHS;
end
resp_LHS = zeros(nLHS,1);
for i = 1:nLHS
    p_LHS = struct();
    for j = 1:nPar
        p_LHS.(sen_param(j).name) = LHSmat(i,j);
    end
    switch resp_ftn
        case 'R'
            resp_LHS(i) = calculateR(p_LHS);
        case {'Cf','Cm','Cmsm'}
            C_end = calculateC(initial,p_LHS,t_end);
            switch resp_ftn
                case 'Cf'
                    resp_LHS(i) = C_end.Cf;
                case 'Cm'
                    resp_LHS(i) = C_end.Cm;
                case 'Cmsm'
                    resp_LHS(i) = C_end.Cmsm;
            end
    end
end
end
[prcc, sign] = PRCC_sen(LHSmat,resp_LHS,idx_influ);
alpha = 0.001;
nInflu = numel(idx_influ);
orientation = 'vertical';
param_names = {sen_param.name};
[~,prccSortedIdx] = sort(abs(prcc),'descend');
prccSortedValue = prcc(prccSortedIdx);
prccSortedLabel = plabel(idx_influ(prccSortedIdx));
prccSortedParamName = param_names(idx_influ(prccSortedIdx));
sig_idx = find(sign < alpha);
for sortIdx = 1:numel(prccSortedIdx)
    fprintf('%s: %f',prccSortedParamName{sortIdx},prccSortedValue(sortIdx));
    if ismember(prccSortedIdx(sortIdx),sig_idx)
        fprintf(' (significant)');
    end
    fprintf('\n');
end
fprintf('\n');
hFig{1}.name = ['sortedPRCC_',resp_ftn];
hFig{1}.figure{1} = figure();
switch orientation
    case 'horizontal'
        plotCount = min(14,numel(prccSortedValue));
        barh(prccSortedValue(1:plotCount))
        set(gca,'YTickLabel',prccSortedLabel(1:plotCount),'YTick',1:plotCount,'FontSize',18)
        title(sprintf('Top %d Significant PRCC',plotCount),'FontSize',18)
        xlabel('PRCC values','FontSize',20)
        ax = gca;
        ax.XAxis.FontSize = 18;
    otherwise
        plotCount = min(8,numel(prccSortedValue));
        b = bar(prccSortedValue(1:plotCount));
        ylabel('Partial rank correlation coefficients (PRCCs)','FontSize',15)
        switch resp_ftn
            case 'Cmsm'
            set(gca,'XTickLabel',prccSortedLabel(1:plotCount),'XTick',1:plotCount,'FontSize',12)
            otherwise
            set(gca,'XTickLabel',prccSortedLabel(1:plotCount),'XTick',1:plotCount,'FontSize',20)
        end
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        ax.TitleFontSizeMultiplier = 1.9;
        ax.YAxis.FontSize = 23;
        set(gcf,'Position',[100,100,1000,560])
        b.BarWidth = 0.6;
end
saveDir = fullfile(projectDir,'outputs','sensitivity',['PRCC_' resp_ftn]);
if ~exist(saveDir,'dir'), mkdir(saveDir); end
T = table(string(param_names(idx_influ(:)))',prcc(:),sign(:), ...
    'VariableNames',{'Parameter','PRCC','PValue'});
writetable(T,fullfile(saveDir,['PRCC_' resp_ftn '.csv']));
for hfIdx=1:numel(hFig)
    exportgraphics(hFig{hfIdx}.figure{1},fullfile(saveDir,[hFig{hfIdx}.name '.png']),'Resolution',300);
    exportgraphics(hFig{hfIdx}.figure{1},fullfile(saveDir,[hFig{hfIdx}.name '.pdf']),'ContentType','vector');
end
end

function repeatPRCC(root,response,seed)
% Independent joint LHS draws with the original OAT selections and ranges.
if strcmp(response,'all'), keys={'R','Cf','Cm','Cmsm'}; else, keys={response}; end
fit=load(fullfile(root,'results','DataFittingVariables.mat'),'initial');
first=load(fullfile(root,'results',['Sensitivity_' keys{1} '.mat']));
par=first.sen_param; n=first.nLHS; X=zeros(n,numel(par));
rng(seed);
for j=1:numel(par)
    X(:,j)=LHS_Call(par(j).range(1),par(j).value,par(j).range(2),0,n,'unif');
end
Y=zeros(n,numel(keys));
for i=1:n
    p=struct();
    for j=1:numel(par), p.(par(j).name)=X(i,j); end
    if any(~strcmp(keys,'R')), cancer=calculateC(fit.initial,p,first.t_end); end
    for j=1:numel(keys)
        if strcmp(keys{j},'R'), Y(i,j)=calculateR(p);
        else, Y(i,j)=cancer.(keys{j}); end
    end
end
T=table();
for j=1:numel(keys)
    ref=load(fullfile(root,'results',['Sensitivity_' keys{j} '.mat']));
    assert(ref.nLHS==n && ref.t_end==first.t_end && ...
        isequal({ref.sen_param.name},{par.name}) && ...
        isequal(vertcat(ref.sen_param.range),vertcat(par.range)), ...
        'Responses must share the same sampling design.');
    [prcc,pval]=PRCC_sen(X,Y(:,j),ref.idx_influ);
    count=numel(prcc); baseRank=zeros(count,1); rank=zeros(count,1);
    [~,order]=sort(abs(ref.prcc),'descend'); baseRank(order)=1:count;
    [~,order]=sort(abs(prcc),'descend'); rank(order)=1:count;
    part=table(repmat(string(keys{j}),count,1),repmat(seed,count,1), ...
        string({par(ref.idx_influ).name})',ref.prcc(:),prcc(:), ...
        prcc(:)-ref.prcc(:),pval(:),baseRank,rank, ...
        'VariableNames',{'Response','Seed','Parameter','BaselinePRCC','PRCC', ...
        'Difference','PValue','BaselineRank','Rank'});
    T=[T;part];
end
saveDir=fullfile(root,'outputs','sensitivity');
if ~exist(saveDir,'dir'),mkdir(saveDir);end
if strcmp(response,'all'), filename=sprintf('PRCC_seed_%d.csv',seed);
else, filename=sprintf('PRCC_%s_seed_%d.csv',response,seed); end
writetable(T,fullfile(saveDir,filename));
fprintf('PRCC seed %d completed using saved OAT selections.\n',seed);
end
