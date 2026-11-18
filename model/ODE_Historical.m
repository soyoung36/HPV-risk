function dydt = ODE_Historical(t,y,p)
% Historical ODE used to generate the 2021 state
p.qH = p.qL;
p.xiH = p.xiL;
yFourGroup = [y(1:32); zeros(6,1)];
dFourGroup = ODE_4groupIncidence(t,yFourGroup,p);
dydt = [dFourGroup(1:32); ...
    dFourGroup(33) + dFourGroup(34); ...
    dFourGroup(35) + dFourGroup(36); ...
    dFourGroup(37) + dFourGroup(38)];
end
