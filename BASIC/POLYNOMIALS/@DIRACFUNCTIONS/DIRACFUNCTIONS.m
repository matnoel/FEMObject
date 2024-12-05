function h = DIRACFUNCTIONS(Q)
% function h = DIRACFUNCTIONS(Xi)
% where Xi is the sample set

h=struct();
param = struct();
param.Q = Q;
domain = 1:Q;
hp = RANDPOLY('DIRACFUNCTIONS',param,domain);
h = class(h,'DIRACFUNCTIONS',hp);
