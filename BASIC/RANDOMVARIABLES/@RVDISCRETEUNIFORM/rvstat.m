function [m,v] = rvstat(u)
param = getparam(u);
[m,v] = unidstat(param.Q);
end
