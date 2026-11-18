function C_end = calculateC(initial,p,t_end)
% Annual new cancer counts
initialValue = [initial.value]';
initialValues = [initialValue; initialValue(end)*0.5];
[~,y] = ode15s(@(t,y) ODE(t,y,p), [0 t_end-1 t_end], initialValues);
C_end = struct();
C_end.Cf = y(end,end-3)-y(end-1,end-3);
C_end.Cm = y(end,end-1)-y(end-1,end-1);
C_end.Cmsm = y(end,end)-y(end-1,end);
end
