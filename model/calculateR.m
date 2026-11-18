function [R,Rmsm]=calculateR(p,v)
if nargin<2
    v.L = (p.xiL+p.qL.*p.mu)./(p.xiL+p.mu);
    v.H = (p.xiH+p.qH.*p.mu)./(p.xiH+p.mu);
    v.M1 = (p.xiM1+p.qM1.*p.mu)./(p.xiM1+p.mu);
    v.M2 = (p.xiM2+p.qM2.*p.mu)./(p.xiM2+p.mu);
end
pL = 1-p.ev.*v.L; pH = 1-p.ev.*v.H;
pM1 = 1-p.ev.*v.M1; pM2 = 1-p.ev.*v.M2;
u1 = p.sigL +p.mu; u2 = p.psi_proL +p.psi_recL +p.mu; u3 = p.g_proL +p.g_recL +p.mu;
u4 = p.sigH +p.mu; u5 = p.psi_proH +p.psi_recH +p.mu; u6 = p.g_proH +p.g_recH +p.mu;
v1 = p.sigM1 +p.mu; v2 = p.psi_proM1 +p.psi_recM1 +p.mu; v3 = p.g_proM1 +p.g_recM1 +p.mu;
v4 = p.sigM2 +p.mu; v5 = p.psi_proM2 +p.psi_recM2 +p.mu; v6 = p.g_proM2 +p.g_recM2 +p.mu;
iL = p.etaL./u1 + p.sigL./(u1.*u2) + p.thetaL.*p.sigL.*p.psi_proL./(u1.*u2.*u3);
iH = p.etaH./u4 + p.sigH./(u4.*u5) + p.thetaH.*p.sigH.*p.psi_proH./(u4.*u5.*u6);
iM1 = p.etaM1./v1 + p.sigM1./(v1.*v2) + p.thetaM1.*p.sigM1.*p.psi_proM1./(v1.*v2.*v3);
iM2 = p.etaM2./v4 + p.sigM2./(v4.*v5) + p.thetaM2.*p.sigM2.*p.psi_proM2./(v4.*v5.*v6);
z13 = p.betMF.*(1-p.aM2).*p.cL.*(pL.*p.piL./p.piM1).*iM1;
z14 = p.betMF.*p.aM2.*p.cL.*(pL.*p.piL./p.piM2).*iM2;
z23 = p.betMF.*p.aH./(1-p.aH).*(1-p.aM2).*p.cL.*(pH.*p.piL./p.piM1).*iM1;
z24 = p.betMF.*p.aH./(1-p.aH).*p.aM2.*p.cL.*(pH.*p.piL./p.piM2).*iM2;
z31 = p.betFM.*(1-p.aM2).*p.cL.*pM1.*iL;
z32 = p.betFM.*p.aH./(1-p.aH).*(1-p.aM2).*p.cL.*(pM1.*p.piL./p.piH).*iH;
z41 = p.betFM.*p.aM2.*p.cL.*pM2.*iL;
z42 = p.betFM.*p.aH./(1-p.aH).*p.aM2.*p.cL.*(pM2.*p.piL./p.piH).*iH;
z44 = p.betMM.*p.cM2M2.*pM2.*iM2;
Rmsm = z44;
all_sizes = [size(z13); size(z14); size(z23); size(z24); ...
             size(z31); size(z32); size(z41); size(z42); size(z44)];
dim = max(all_sizes);
R = zeros(dim);
for idx1 = 1:dim(1)
    for idx2 = 1:dim(2)
        Z = zeros(4);
        Z(1,3)=z13(min(idx1,end),min(idx2,end)); Z(1,4)=z14(min(idx1,end),min(idx2,end));
        Z(2,3)=z23(min(idx1,end),min(idx2,end)); Z(2,4)=z24(min(idx1,end),min(idx2,end));
        Z(3,1)=z31(min(idx1,end),min(idx2,end)); Z(3,2)=z32(min(idx1,end),min(idx2,end));
        Z(4,1)=z41(min(idx1,end),min(idx2,end)); Z(4,2)=z42(min(idx1,end),min(idx2,end));
        Z(4,4)=z44(min(idx1,end),min(idx2,end));
        R(idx1,idx2) = max(eig(Z));
    end
end
end
