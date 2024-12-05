function [L,P]=orthogonalize(L)


if all(size(L)>1)
    error('doit etre un vecteur de variable aleatoires')
end

if nargin==1
    tol=100*eps;
end

D = expectmtimes(L,L');

s = eig(full(D));
rang = sum(sqrt(abs(s/max(s)))>tol);

%keyboard
opts.disp=0;

if rang<numel(L)
    [V,S] = eigs(D,rang,'LM',opts);
else
    [V,S] = eig(full(D));    
end


P = V';
L = P*L;
