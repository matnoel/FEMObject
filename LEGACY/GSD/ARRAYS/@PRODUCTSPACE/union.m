function [c,Ia,Ib] = union(a,b)
% function [c,Ia,Ib] = union(a,b)

dima = getdim(a);
dimb = getdim(b);
[dim,Ia,Ib] = union(dima,dimb);
c = PRODUCTSPACE(dim);




