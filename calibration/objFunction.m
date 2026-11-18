function res = objFunction(x_fit, data, initial, param)
    fitIdx = find([param.fit_mask]);
    for i = 1:length(x_fit)
        param(fitIdx(i)).value = x_fit(i);
    end
    rho_pro = param(strcmp({param.name}, 'rho_pro')).value;
    rho_rec = param(strcmp({param.name}, 'rho_rec')).value;
    map_pro = { 'psi_proL', 'psi_proL', 'g_proL', 'g_proL';
                'psi_proM1','psi_proM2','g_proM1','g_proM2' };
    map_rec = { 'psi_recL','psi_recL','g_recL','g_recL';
                'psi_recM1','psi_recM2','g_recM1','g_recM2' };
    for k = 1:size(map_pro,2)
        f_idx = strcmp({param.name}, map_pro{1,k});
        m_idx = strcmp({param.name}, map_pro{2,k});
        param(m_idx).value = rho_pro * param(f_idx).value;
    end
    for k = 1:size(map_rec,2)
        f_idx = strcmp({param.name}, map_rec{1,k});
        m_idx = strcmp({param.name}, map_rec{2,k});
        param(m_idx).value = rho_rec * param(f_idx).value;
    end
    for cancerIdx = find(ismember({param.name}, {'wM1','wM2'}))
        param(cancerIdx).value = 0.45;
    end
    map_f = { 'qL','zetL','etaL','thetaL','xiL','sigL','psi_proL','psi_recL','g_proL','g_recL','wL','delL'; ...
              'qH','zetH','etaH','thetaH','xiH','sigH','psi_proH','psi_recH','g_proH','g_recH','wH','delH' };
    for k = 1:size(map_f,2)
        l_idx = strcmp({param.name}, map_f{1,k});
        h_idx = strcmp({param.name}, map_f{2,k});
        param(h_idx).value = param(l_idx).value;
    end
    tspan = data.time.all-data.time.all(1);
    sim = ODEsim_cancer(tspan, initial, param);
    resCCf = (data.annual.f' - gradient(sim.CCf))/ max(data.annual.f);
    resCDf = (data.annual.f_death' - gradient(sim.CDf))/ max(data.annual.f_death);
    resCCm = (data.annual.m' - gradient(sim.CCm))/ max(data.annual.m);
    res = [resCCf; resCDf; resCCm];
end
