function [c,I] = setdiff(a,b)
dima = getdim(a);
dimb = getdim(b);
[dim,I] = setdiff(dima,dimb);
c = PRODUCTSPACE(dim);





