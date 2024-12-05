function [s,I] = ismember(a,b)
%function [s,I] = ismember(a,b)

dima = getdim(a);
dimb = getdim(b);

[s,I] = ismember(dima,dimb);


