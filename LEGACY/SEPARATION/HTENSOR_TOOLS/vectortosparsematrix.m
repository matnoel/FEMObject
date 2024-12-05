function v = vectortosparsematrix(v,P,n,varargin)
% function H = vectortosparsematrix(v,P,n)
% function H = vectortosparsematrix(v,P,n,sym_dim)


if nargin == 3
    dim_sym = [];
else
    dim_sym = varargin{1};
end

dim = ndims(v);
d2i = v.dim2ind;
r = rank(v);

for d = dim_sym
    i = d2i(d);
    nn = prod(n(:,d));

    [I,J] = ind2sub(n(:,d),P{d});
    ndiag = ~(I==J);
    vUi = [v.U{i}; v.U{i}(ndiag,:)];

    I1 = [I;J(ndiag)];
    J1 = [J;I(ndiag)];
    I = sub2ind(n(:,d),I1,J1);

    I = repmat(I,r(i),1);
    J = ones(size(I1));
    J = J*(1:r(i));
    v.U{i} = sparse(I,J(:),vUi,nn,r(i));
end

for d = setdiff(1:dim,dim_sym)
    i = d2i(d);
    nn = prod(n(:,d));
    I = repmat(P{d},r(i),1);
    J = ones(size(P{d}));
    J = J*(1:r(i));
    v.U{i} = sparse(I,J(:),v.U{i},nn,r(i));
end

end
