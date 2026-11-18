function data=createData(dataFile)
% Read processed annual cancer data (1999-2022)
raw = readmatrix(dataFile);
data.time.all = raw(1,2:end);
data.annual.f = raw(2,2:end);
data.annual.f_death = raw(3,2:end);
data.annual.m = raw(4,2:end);
data.cumul.f = cumsum(data.annual.f);
data.cumul.f_death = cumsum(data.annual.f_death);
data.cumul.m = cumsum(data.annual.m);
end
