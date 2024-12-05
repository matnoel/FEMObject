function A = normalizefactors(A,dim)
% function A = normalizefactors(A,dim)

if nargin==1
    dim=1:A.dim;
end

ndim = numel(dim);

for i=1:A.m
    if A.alpha(i) < 0
        A.F{i,1} = -A.F{i,1};
        A.alpha(i) = -A.alpha(i);
    end
    salpha = A.alpha(i)^(1/ndim);
    A.F(i,dim) = cellfun(@(x) salpha*x,A.F(i,dim),'uniformoutput',0);
    A.alpha(i)=1;
end
