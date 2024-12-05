function P = getpattern(A,varargin)
% function P = getpattern(A)
% function P = getpattern(A,dim_sym)
% dim_sym : symmetric dimensions

dim=A.dim;
m=A.m;

P=cell(1,dim);

if nargin == 2
    dim_sym=varargin{1};
else
    dim_sym=[];
end

for d = dim_sym
    P{d}=[];
    for i=1:m
        P{d} = union(P{d},find(triu(A.F{i,d})));
    end
end

for d = setdiff(1:dim,dim_sym)
    P{d}=[];
    for i=1:m
        P{d} = union(P{d},find(A.F{i,d}));
    end
end

end
