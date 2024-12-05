function [c,Ia,Ib] = intersect(a,b)
% function [c,Ia,Ib] = intersect(a,b)

dima = getdim(a);
dimb = getdim(b);
[dim,Ia,Ib] = intersect(dima,dimb);
c = PRODUCTSPACE(dim);





