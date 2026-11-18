function sim = ODEsim_cancer(tspan, initial, param)
    for i = 1:numel(param)
        p.(param(i).name) = param(i).value;
    end
    p.cM2M2 = p.cM2*p.rM2;
    [~, ySim] = solvePiecewiseODE(tspan, initial, p);
    sim.CCf = ySim(:, end-2);
    sim.CDf = ySim(:, end-1);
    sim.CCm = ySim(:, end);
end
