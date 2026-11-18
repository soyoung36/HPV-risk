function dydt=ODE_4groupIncidence(t,y,p)
% Projection ODE with group-specific cumulative case/death trackers
Sl  = y(1); Vl  = y(2); El  = y(3); Il  = y(4); Pl  = y(5); Cl  = y(6); Rcl = y(7); Rl  = y(8);
Sh  = y(9); Vh  = y(10); Eh  = y(11); Ih  = y(12); Ph  = y(13); Ch  = y(14); Rch = y(15); Rh  = y(16);
Sm1  = y(17); Vm1  = y(18); Em1  = y(19); Im1  = y(20); Pm1  = y(21); Cm1  = y(22); Rcm1 = y(23); Rm1  = y(24);
Sm2  = y(25); Vm2  = y(26); Em2  = y(27); Im2  = y(28); Pm2  = y(29); Cm2  = y(30); Rcm2 = y(31); Rm2  = y(32);
Nl = sum(y(1:8)); Nh = sum(y(9:16)); Nm1 = sum(y(17:24)); Nm2 = sum(y(25:32));
map_pro = { 'psi_proL', 'psi_proL', 'g_proL', 'g_proL';
            'psi_proM1','psi_proM2','g_proM1','g_proM2' };
map_rec = { 'psi_recL','psi_recL','g_recL','g_recL';
            'psi_recM1','psi_recM2','g_recM1','g_recM2' };
for k = 1:size(map_pro,2)
    p.(map_pro{2,k}) = p.rho_pro * p.(map_pro{1,k});
end
for k = 1:size(map_rec,2)
    p.(map_rec{2,k}) = p.rho_rec * p.(map_rec{1,k});
end
p.wM1 = 0.45; p.wM2 = 0.45;
p.zetH = p.zetL; p.etaH = p.etaL; p.thetaH = p.thetaL; p.sigH = p.sigL;
p.psi_proH = p.psi_proL; p.psi_recH = p.psi_recL; p.g_proH = p.g_proL; p.g_recH = p.g_recL; p.wH = p.wL; p.delH = p.delL;
jl = p.etaL*El +Il +p.thetaL*Pl;
jh = p.etaH*Eh +Ih +p.thetaH*Ph;
jm1 = p.etaM1*Em1 +Im1 +p.thetaM1*Pm1;
jm2 = p.etaM2*Em2 +Im2 +p.thetaM2*Pm2;
lamL = p.betMF*p.cL* ((1-p.aM2)*jm1/Nm1 +p.aM2*jm2/Nm2);
lamH = p.betMF*p.cL* (p.aH/(1-p.aH)*(1-p.aM2)*Nl/Nh*jm1/Nm1 +p.aH/(1-p.aH)*p.aM2*Nl/Nh*jm2/Nm2);
lamM1 = p.betFM*p.cL*(1-p.aM2)* (jl/Nm1 +p.aH/(1-p.aH)*Nl/Nm1*jh/Nh);
lamM2 = p.betFM*p.cL*p.aM2* (jl/Nm2 +p.aH/(1-p.aH)*Nl/Nm2*jh/Nh) +p.betMM*p.cM2M2*jm2/Nm2;
dSl = (1-p.qL)*p.piL +p.zetL*Rl -(lamL +p.xiL +p.mu)*Sl;
dVl = p.qL*p.piL +p.xiL*Sl -(1-p.ev)*lamL*Vl -p.mu*Vl;
dEl = lamL*Sl +(1-p.ev)*lamL*Vl -(p.sigL +p.mu)*El;
dIl = p.sigL*El -(p.psi_proL +p.psi_recL +p.mu)*Il;
dPl = p.psi_proL*Il -(p.g_proL +p.g_recL +p.mu)*Pl;
dCl = p.g_proL*Pl -(p.wL +p.delL +p.mu)*Cl;
dRcl = p.wL*Cl -p.mu*Rcl;
dRl = p.psi_recL*Il +p.g_recL*Pl -(p.zetL +p.mu)*Rl;
dSh = (1-p.qH)*p.piH +p.zetH*Rh -(lamH +p.xiH +p.mu)*Sh;
dVh = p.qH*p.piH +p.xiH*Sh -(1-p.ev)*lamH*Vh -p.mu*Vh;
dEh = lamH*Sh +(1-p.ev)*lamH*Vh -(p.sigH +p.mu)*Eh;
dIh = p.sigH*Eh -(p.psi_proH +p.psi_recH +p.mu)*Ih;
dPh = p.psi_proH*Ih -(p.g_proH +p.g_recH +p.mu)*Ph;
dCh = p.g_proH*Ph -(p.wH +p.delH +p.mu)*Ch;
dRch = p.wH*Ch -p.mu*Rch;
dRh = p.psi_recH*Ih +p.g_recH*Ph -(p.zetH +p.mu)*Rh;
dSm1 = (1-p.qM1)*p.piM1 +p.zetM1*Rm1 -(lamM1 +p.xiM1 +p.mu)*Sm1;
dVm1 = p.qM1*p.piM1 +p.xiM1*Sm1 -(1-p.ev)*lamM1*Vm1 -p.mu*Vm1;
dEm1 = lamM1*Sm1 +(1-p.ev)*lamM1*Vm1 -(p.sigM1 +p.mu)*Em1;
dIm1 = p.sigM1*Em1 -(p.psi_proM1 +p.psi_recM1 +p.mu)*Im1;
dPm1 = p.psi_proM1*Im1 -(p.g_proM1 +p.g_recM1 +p.mu)*Pm1;
dCm1 = p.g_proM1*Pm1 -(p.wM1 +p.delM1 +p.mu)*Cm1;
dRcm1 = p.wM1*Cm1 -p.mu*Rcm1;
dRm1 = p.psi_recM1*Im1 +p.g_recM1*Pm1 -(p.zetM1 +p.mu)*Rm1;
dSm2 = (1-p.qM2)*p.piM2 +p.zetM2*Rm2 -(lamM2 +p.xiM2 +p.mu)*Sm2;
dVm2 = p.qM2*p.piM2 +p.xiM2*Sm2 -(1-p.ev)*lamM2*Vm2 -p.mu*Vm2;
dEm2 = lamM2*Sm2 +(1-p.ev)*lamM2*Vm2 -(p.sigM2 +p.mu)*Em2;
dIm2 = p.sigM2*Em2 -(p.psi_proM2 +p.psi_recM2 +p.mu)*Im2;
dPm2 = p.psi_proM2*Im2 -(p.g_proM2 +p.g_recM2 +p.mu)*Pm2;
dCm2 = p.g_proM2*Pm2 -(p.wM2 +p.delM2 +p.mu)*Cm2;
dRcm2 = p.wM2*Cm2 -p.mu*Rcm2;
dRm2 = p.psi_recM2*Im2 +p.g_recM2*Pm2 -(p.zetM2 +p.mu)*Rm2;
dCCl = p.g_proL*Pl;
dCCh = p.g_proH*Ph;
dCDl = p.delL*Cl;
dCDh = p.delH*Ch;
dCCm1 = p.g_proM1*Pm1;
dCCm2 = p.g_proM2*Pm2;
dydt = [dSl; dVl; dEl; dIl; dPl; dCl; dRcl; dRl; ...
        dSh; dVh; dEh; dIh; dPh; dCh; dRch; dRh; ...
        dSm1; dVm1; dEm1; dIm1; dPm1; dCm1; dRcm1; dRm1; ...
        dSm2; dVm2; dEm2; dIm2; dPm2; dCm2; dRcm2; dRm2; ...
        dCCl; dCCh; dCDl; dCDh; dCCm1; dCCm2];
end
