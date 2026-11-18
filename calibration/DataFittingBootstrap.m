function DataFittingBootstrap(fresh,nRequested)
% Replot saved bootstrap draws, or refit residual bootstrap draws when fresh=true.
% Annual bands use percentiles of annual projections from each bootstrap draw.
if nargin<1, fresh=false; end
if nargin<2, nRequested=1000; end
% Load the baseline fit and set the bootstrap options.
set(0, 'defaultaxesfontsize',18, 'defaulttextinterpreter','latex')
projectDir = setup_paths();
rng('default')
resultFile = fullfile(projectDir, 'results', 'DataFittingVariables.mat');
load(resultFile)
saveDir=fullfile(projectDir,'outputs','calibration');
if ~exist(saveDir,'dir'), mkdir(saveDir); end
fields = {'CCf','CDf','CCm'};

% Resample centered residuals and refit, or load the saved estimates.
if fresh
    nBoot = nRequested;
    nFitData = numel(data.time.all);
    bootResiduals = [];
    for i = 1:3
        residual_separate.(fields{i}) = residual((nFitData*(i-1)+1):(nFitData*i));
        residual_separate.(fields{i}) = residual_separate.(fields{i}) ...
            - mean(residual_separate.(fields{i}));
        [~, bootIdx] = bootstrp(nBoot, [], residual_separate.(fields{i}));
        bootResiduals_separate.(fields{i}) = residual_separate.(fields{i})(bootIdx);
        bootResiduals = [bootResiduals; bootResiduals_separate.(fields{i})];
    end
    fitIdx = find([final_param.fit_mask]);
    x0 = [final_param(fitIdx).value];
    lb = arrayfun(@(x) x.range(1), final_param(fitIdx));
    ub = arrayfun(@(x) x.range(2), final_param(fitIdx));
    boot_x = zeros(nBoot,length(fitIdx));
    boot_resnorm = zeros(nBoot,1);
    boot_exitflag = zeros(nBoot,1);
    options = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-12, 'MaxFunEvals',1e3);
    tic;
    for i = 1:nBoot
        fprintf('Running %d th Bootstrap Iterations...\n', i);
        bootData = data;
        bootData.annual.f = gradient(sim.CCf(1:nFitData))' + bootResiduals(1:nFitData,i)' * max(data.annual.f);
        bootData.annual.f_death = gradient(sim.CDf(1:nFitData))' + bootResiduals(nFitData+1:2*nFitData,i)' * max(data.annual.f_death);
        bootData.annual.m = gradient(sim.CCm(1:nFitData))' + bootResiduals(2*nFitData+1:end,i)'* max(data.annual.m);
        [boot_x(i,:), boot_resnorm(i), ~, boot_exitflag(i)] = lsqnonlin(...
            @(x) objFunction(x, bootData, initial, param), x0, lb, ub, options);
    end
    toc;
    CI_95 = prctile(boot_x, [2.5, 50, 97.5]);
    fprintf('\n--- Bootstrap Result (95%% intervals) ---\n');
    for i = 1:length(fitIdx)
        fprintf('%s: Fitted=%.4f, Bootstrap median=%.4f, 95%% CI=[%.4f, %.4f]\n', ...
            final_param(fitIdx(i)).name, final_param(fitIdx(i)).value, ...
            CI_95(2,i), CI_95(1,i), CI_95(3,i));
    end
    bootstrapResultFile = fullfile(projectDir, 'outputs', 'calibration', 'DataFittingBootstrapResult.mat');
    save(bootstrapResultFile, 'boot_x', 'CI_95', 'boot_resnorm', 'boot_exitflag')
else
    load(fullfile(projectDir,'results','DataFittingBootstrapResult.mat'),'boot_x','CI_95');
end

% Propagate the parameter estimates through the model.
set(0, 'defaultaxesfontsize',18, 'defaulttextinterpreter','latex')
nBoot = size(boot_x,1);
fitIdx = find([final_param.fit_mask]);
tspan = data.time.all-data.time.all(1);
boot_sims = zeros(length(tspan), 3, nBoot);
boot_R = zeros(1, nBoot);
for i = 1:nBoot
    temp_param = final_param;
    for j = 1:length(fitIdx)
        temp_param(fitIdx(j)).value = boot_x(i, j);
    end
    rho_pro = temp_param(strcmp({temp_param.name}, 'rho_pro')).value;
    rho_rec = temp_param(strcmp({temp_param.name}, 'rho_rec')).value;
    map_pro = { 'psi_proL', 'psi_proL', 'g_proL', 'g_proL';
                'psi_proM1','psi_proM2','g_proM1','g_proM2' };
    map_rec = { 'psi_recL','psi_recL','g_recL','g_recL';
                'psi_recM1','psi_recM2','g_recM1','g_recM2' };
    for k = 1:size(map_pro,2)
        f_idx = strcmp({temp_param.name}, map_pro{1,k});
        m_idx = strcmp({temp_param.name}, map_pro{2,k});
        temp_param(m_idx).value = rho_pro * temp_param(f_idx).value;
    end
    for k = 1:size(map_rec,2)
        f_idx = strcmp({temp_param.name}, map_rec{1,k});
        m_idx = strcmp({temp_param.name}, map_rec{2,k});
        temp_param(m_idx).value = rho_rec * temp_param(f_idx).value;
    end
    for cancerIdx = find(ismember({temp_param.name}, {'wM1','wM2'}))
        temp_param(cancerIdx).value = 0.45;
    end
    map_f = { 'qL','zetL','etaL','thetaL','xiL','sigL','psi_proL','psi_recL','g_proL','g_recL','wL','delL'; ...
              'qH','zetH','etaH','thetaH','xiH','sigH','psi_proH','psi_recH','g_proH','g_recH','wH','delH' };
    for k = 1:size(map_f,2)
        l_idx = strcmp({temp_param.name}, map_f{1,k});
        h_idx = strcmp({temp_param.name}, map_f{2,k});
        temp_param(h_idx).value = temp_param(l_idx).value;
    end
    sim_boot = ODEsim_cancer(tspan, initial, temp_param);
    boot_sims(:,:,i) = [sim_boot.CCf, sim_boot.CDf, sim_boot.CCm];
    for j = 1:numel(temp_param)
        temp_p.(temp_param(j).name) = temp_param(j).value;
    end
    temp_p.cM2M2 = temp_p.cM2*temp_p.rM2;
    boot_R(i) = calculateR(temp_p);
end
lower_bound = prctile(boot_sims, 2.5, 3);
upper_bound = prctile(boot_sims, 97.5, 3);
annual_boot = zeros(size(boot_sims));
for j=1:3
    for k=1:nBoot
        annual_boot(:,j,k) = gradient(boot_sims(:,j,k));
    end
end
annual_lower = prctile(annual_boot,2.5,3);
annual_upper = prctile(annual_boot,97.5,3);
CI_R = prctile(boot_R, [2.5, 50, 97.5]);
fprintf('Reproduction number: Bootstrap median=%.4f, 95%% CI=[%.4f, %.4f]\n', CI_R(2), CI_R(1), CI_R(3));

% Plot the bootstrap bands and save
titles = {'(a) Annual Cervical Cancer Cases', '(b) Annual Cervical Cancer Deaths', '(c) Annual Male Cancer Cases',...
    '(d) Cumulative Cervical Cancer Cases', '(e) Cumulative Cervical Cancer Deaths', '(f) Cumulative Male Cancer Cases'};
data_fields = {'f', 'f_death', 'm'};
figure('Position', [100, 200, 1800, 1000]);
for i = 1:3
    subplot(2, 3, i);
    hold on;
    fill([data.time.all'; flipud(data.time.all')],[annual_lower(:,i);flipud(annual_upper(:,i))],[0 0.4471 0.7412],'EdgeColor', 'none', 'FaceAlpha', 0.45);
    hold on
    p1 = plot(data.time.all, gradient(sim.(fields{i})), 'b', 'LineWidth', 3);
    data_i = [data.annual.(data_fields{i})];
    p4 = scatter(data.time.all, data_i, 80, 'r', 'filled');
    title(titles{i});
    xlabel('Time (year)');
    ylabel('Number of Cases');
    grid on;
    if i == 1
        legend([p1,p4], {'Model Fit', 'Data'}, 'Location', 'best');
    end
    if i == 2
        ylabel('Number of Deaths')
    end
    ylim([0 inf])
end
for i = 1:3
    subplot(2, 3, i+3);
    hold on;
    fill([data.time.all'; flipud(data.time.all')],[lower_bound(:,i);flipud(upper_bound(:,i))],[0 0.4471 0.7412],'EdgeColor', 'none', 'FaceAlpha', 0.45);
    hold on
    p1 = plot(data.time.all, sim.(fields{i}), 'b', 'LineWidth', 3);
    data_i = [data.cumul.(data_fields{i})];
    p4 = scatter(data.time.all, data_i, 80, 'r', 'filled');
    title(titles{i+3});
    xlabel('Time (year)');
    ylabel('Number of Cases');
    grid on;
    ax = gca;
    ax.YAxis.Exponent = 0;
    if i == 1
        legend([p1,p4], {'Model Projection', 'Data'}, 'Location', 'southeast');
    end
    if i == 2
        ylabel('Number of Deaths')
    end
end
saveDir = fullfile(projectDir, 'outputs', 'calibration');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
exportgraphics(gcf, fullfile(saveDir, 'Fitting_with_CI.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(saveDir, 'Fitting_with_CI.pdf'), 'ContentType', 'vector');
end
