function s = eq(a,b)
% function s = eq(a,b)

dima = getdim(a);
dimb = getdim(b);
if length(dima)==length(dimb) && all(dima==dimb)
    s = true;
else
    s = false;
end

