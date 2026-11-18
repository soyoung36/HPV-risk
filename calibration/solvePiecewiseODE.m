function [tSim, ySim] = solvePiecewiseODE(tspan, initial, p)
% Integrate across the vaccination switch at t=8 (2007)
tspan = tspan(:)';
switchTime = 2007 - 1999;
y0 = [initial.value]';
pPre = p;
pPre.qL = 0;
pPre.qH = 0;
pPre.xiL = 0;
pPre.xiH = 0;
pPost = p;
pPost.qL = 0.05;
pPost.qH = 0.05;
pPost.xiL = 0.01;
pPost.xiH = 0.01;
if tspan(1) ~= 0
    error('solvePiecewiseODE:InvalidStartTime', ...
        'tspan must begin at t = 0, corresponding to 1999.');
end
if tspan(end) <= switchTime
    [tSim, ySim] = ode45(@(t,y) ODE_Historical(t,y,pPre), tspan, y0);
    return
end
preTimes = unique([tspan(tspan < switchTime), switchTime]);
[tPre, yPre] = ode45(@(t,y) ODE_Historical(t,y,pPre), preTimes, y0);
postTimes = unique([switchTime, tspan(tspan > switchTime)]);
[tPost, yPost] = ode45(@(t,y) ODE_Historical(t,y,pPost), postTimes, yPre(end,:)');
tCombined = [tPre; tPost(2:end)];
yCombined = [yPre; yPost(2:end,:)];
[isPresent, outputIdx] = ismember(tspan, tCombined);
if ~all(isPresent)
    error('solvePiecewiseODE:MissingOutputTime', ...
        'The piecewise solution does not contain every requested time.');
end
tSim = tCombined(outputIdx);
ySim = yCombined(outputIdx,:);
end
