function P = getpattern(H,varargin)
% function P = getpattern(H)
% function P = getpattern(H,dim_sym)
% dim_sym : symmetric dimensions

dim = ndims(H);
d2i = H.dim2ind;
r = rank(H);

P = cell(1,dim);

if nargin == 2
    dim_sym = varargin{1};
else
    dim_sym = [];
end

for d = dim_sym
    P{d} = [];
    i = d2i(d);
    n = sqrt(size(H.U{i},1));
    for m = 1:r(i)
        HUm = H.U{i}(:,m);
        HUm = reshape(HUm,n,n);
        HUm = triu(HUm);
        P{d} = union(P{d},find(HUm));
    end
end

for d = setdiff(1:dim,dim_sym)
    P{d} = [];
    i = d2i(d);
    for m = 1:r(i)
        P{d} = union(P{d},find(H.U{i}(:,m)));
    end
end

end
