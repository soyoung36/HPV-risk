function initial=createInitial(initialData)
% Initial states and cumulative trackers
raw_initial = {
    'Sl', 11732500;
    'Vl', 0;
    'El', 83600;
    'Il', 551000;
    'Pl', 95000;
    'Cl', 3591.2;
    'Rcl', 0;
    'Rl', 570000;
    'Sh', 617500;
    'Vh', 0;
    'Eh', 4400;
    'Ih', 29000;
    'Ph', 5000;
    'Ch', 897.8;
    'Rch', 0;
    'Rh', 30000;
    'Sm1', 11640000;
    'Vm1', 0;
    'Em1', 116400;
    'Im1', 194000;
    'Pm1', 11640;
    'Cm1', 111.965;
    'Rcm1', 0;
    'Rm1', 776000;
    'Sm2', 260000;
    'Vm2', 0;
    'Em2', 3600;
    'Im2', 6000;
    'Pm2', 360;
    'Cm2', 111.965;
    'Rcm2', 0;
    'Rm2', 124000;
    'CCf', initialData.Cf0;
    'CDf', initialData.Df0;
    'CCm', initialData.Cm0;
};
nState = size(raw_initial, 1);
initial = struct('state', {}, 'value', {});
for i = 1:nState
    initial(i).state  = raw_initial{i, 1};
    initial(i).value = raw_initial{i, 2};
end
end
