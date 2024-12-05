function [x,n] = normalizephi(x,dim)
if nargin==1
    dim = 1:getnbdim(x);
end
n=1;
for i=dim
    nt = norm(x.phi{i});
    n = n*nt;
    x.phi{i} = x.phi{i}/nt;
end

