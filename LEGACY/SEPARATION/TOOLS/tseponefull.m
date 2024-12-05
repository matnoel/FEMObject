function u = tseponefull(n)
% function u = tseponefull(n)

dim = length(n);
m = n;
if dim==1
    alpha=tensor(ones(m,1),m);
else
    alpha = tensor(ones(m),m);
end
F = cell(dim,1);
for i=1:dim
    F{i}=speye(m(i));
end

u = splitvectors(TSEPMATRIX(F,alpha));
