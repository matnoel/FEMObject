function [a,xi] = random(pcr,n)
% function [a,xi] = random(pcr,n)

if nargin==1
    n = 1;
end

L = pcr.D * pcr.L;
[L,xi] = random(L,1,prod(n));
if length(n)==1
    n = [n,1];
end
if prod(n)>1
    L = MULTIMATRIX(L,[pcr.m,1],n);
end

if issparse(pcr.V)
    pcr.V = sparse(pcr.V);
    L = sparse(L);
end

if iscell(pcr.V)
    a = pcr.V{1}*L(1);
    for k=2:pcr.m
        a = a + pcr.V{k}*L(k);
    end
else
    a = double(pcr.V)*L;
end

if all(size(a)==1)
    a = reshape(double(a),n);
else
    if prod(n)==1
        a = reshape(double(a),size(pcr));
    else
        a = reshape(a,size(pcr));
    end
end
