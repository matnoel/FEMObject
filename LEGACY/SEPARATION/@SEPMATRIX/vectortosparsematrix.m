function v = vectortosparsematrix(v,P,n,varargin)
% function S = vectortosparsematrix(v,P,n)
% function S = vectortosparsematrix(v,P,n,sym_dim)

if nargin == 3
    dim_sym = [];
else
    dim_sym = varargin{1};
end

for d = dim_sym
    [I,J] = ind2sub(n(:,d),P{d});
    ndiag = ~(I==J);
    I1 = [I;J(ndiag)];
    J1 = [J;I(ndiag)];
    for m = 1:v.m
        F = [v.F{m,d}; v.F{m,d}(ndiag)];
        v.F{m,d} = sparse(I1,J1,F,...
            n(1,d),n(2,d));
    end
end

for d = setdiff(1:v.dim,dim_sym)
    [I,J] = ind2sub(n(:,d),P{d});
    for m = 1:v.m
        v.F{m,d} = sparse(I,J,v.F{m,d},...
            n(1,d),n(2,d));
    end
end
