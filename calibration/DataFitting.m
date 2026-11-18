function DataFitting()
% Load data, initial states, and parameters
close all
set(0, 'defaultaxesfontsize',20, 'defaulttextinterpreter','latex')
projectDir = setup_paths();
dataFile = fullfile(projectDir, 'data', 'proc_KoreaHPVcancer_1999-2022.csv');
data = createData(dataFile);
plotData(data)
saveDir = fullfile(projectDir, 'outputs', 'calibration');
if ~exist(saveDir, 'dir'), mkdir(saveDir); end
exportgraphics(gcf, fullfile(saveDir, 'Data.png'), 'Resolution', 300);
initialData.Cf0 = data.annual.f(1);
initialData.Df0 = data.annual.f_death(1);
initialData.Cm0 = data.annual.m(1);
initial = createInitial(initialData);
param = createParam();
target_params = {'rho_pro','rho_rec'};
for i = 1:length(param)
    if ismember(param(i).name, target_params)
        param(i).fit_mask = true;
    end
end
fitIdx = find([param.fit_mask]);
x0 = [param(fitIdx).value];
lb = arrayfun(@(x) x.range(1), param(fitIdx));
ub = arrayfun(@(x) x.range(2), param(fitIdx));

% Fit parameters
close all
options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12, 'MaxFunEvals',1e3);
[best_x, resnorm, residual, exitflag] = lsqnonlin(...
    @(x) objFunction(x, data, initial, param), x0, lb, ub, options);

% Apply the fitted ratios and simulate the model
final_param = param;
for i = 1:length(best_x)
    final_param(fitIdx(i)).value = best_x(i);
    fprintf('Optimized %s: %.4f\n', final_param(fitIdx(i)).name, best_x(i));
end
rho_pro = final_param(strcmp({final_param.name}, 'rho_pro')).value;
rho_rec = final_param(strcmp({final_param.name}, 'rho_rec')).value;
map_pro = { 'psi_proL', 'psi_proL', 'g_proL', 'g_proL';
            'psi_proM1','psi_proM2','g_proM1','g_proM2' };
map_rec = { 'psi_recL','psi_recL','g_recL','g_recL';
            'psi_recM1','psi_recM2','g_recM1','g_recM2' };
for k = 1:size(map_pro,2)
    f_idx = strcmp({final_param.name}, map_pro{1,k});
    m_idx = strcmp({final_param.name}, map_pro{2,k});
    final_param(m_idx).value = rho_pro * final_param(f_idx).value;
end
for k = 1:size(map_rec,2)
    f_idx = strcmp({final_param.name}, map_rec{1,k});
    m_idx = strcmp({final_param.name}, map_rec{2,k});
    final_param(m_idx).value = rho_rec * final_param(f_idx).value;
end
for cancerIdx = find(ismember({final_param.name}, {'wM1','wM2'}))
    final_param(cancerIdx).value = 0.45;
end
map_f = { 'qL','zetL','etaL','thetaL','xiL','sigL','psi_proL','psi_recL','g_proL','g_recL','wL','delL'; ...
          'qH','zetH','etaH','thetaH','xiH','sigH','psi_proH','psi_recH','g_proH','g_recH','wH','delH' };
for k = 1:size(map_f,2)
    l_idx = strcmp({final_param.name}, map_f{1,k});
    h_idx = strcmp({final_param.name}, map_f{2,k});
    final_param(h_idx).value = final_param(l_idx).value;
end
residual = objFunction(best_x, data, initial, final_param);
for i = 1:numel(final_param)
    final_p.(final_param(i).name) = final_param(i).value;
end
final_p.cM2M2 = final_p.cM2*final_p.rM2;
Rv = calculateR(final_p);
fprintf('Reproduction number: %g\n',Rv);
tspan = (data.time.all-data.time.all(1));
sim = ODEsim_cancer(tspan,initial,final_param);

% Plot the results and save
titles = {'Annual Female Cases', 'Annual Female Deaths', 'Annual Male Cases',...
    'Cumulative Female Cases', 'Cumulative Female Deaths', 'Cumulative Male Cases'};
sim_fields = {'CCf', 'CDf', 'CCm'};
data_fields = {'f', 'f_death', 'm'};
fit_metrics = struct();
for i = 1:3
    observed = data.annual.(data_fields{i});
    predicted = gradient(sim.(sim_fields{i}));
    observed = observed(:);
    predicted = predicted(:);
    metric.n = numel(observed);
    metric.mean_observed = mean(observed);
    metric.rmse = sqrt(mean((observed - predicted).^2));
    metric.nrmse = metric.rmse / metric.mean_observed;
    fit_metrics.(data_fields{i}) = metric;
    fprintf('%s: RMSE = %.4f; mean-normalized RMSE = %.6f (%.2f%%), n = %d\n', ...
        titles{i}, metric.rmse, metric.nrmse, 100*metric.nrmse, metric.n);
end
resultFile = fullfile(projectDir, 'outputs', 'calibration', 'DataFittingVariables.mat');
save(resultFile, 'data', 'initial', 'param', 'final_param', 'residual', 'sim', 'fit_metrics', 'exitflag', 'resnorm')
end
