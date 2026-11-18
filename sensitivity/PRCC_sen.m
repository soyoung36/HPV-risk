function [prcc,sign] = PRCC_sen(LHSmat,resp_LHS,idx_influ)
% Adapted from the Marino/Renardy LHS-PRCC code. See THIRD_PARTY_NOTICES.md
Y = resp_LHS(:);
nInflu = numel(idx_influ);
prcc = zeros(1,nInflu);
sign = zeros(1,nInflu);
for j = 1:nInflu
    i = idx_influ(j);
    Z = LHSmat;
    Z(:,i) = [];
    [rho,pval] = partialcorr([LHSmat(:,i),Y],Z,'type','Spearman');
    prcc(j) = rho(1,2);
    sign(j) = pval(1,2);
end
end
