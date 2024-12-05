function r=uniqueV(r,tol)
% function r=uniqueV(r,tol)
if size(r,2)>1
 warning('les variables deterministes sont des matrices : elles sont selectionnees comme des vecteurs')
end

if nargin==1
    tol=100*eps;
end

V = double(r.V);
s = svds(V,size(V,2));
rang = sum(s/max(s)>tol);


if rang<r.m

[V,S,U] = svds(V,rang);

r.m = rang;

r.V = MULTIMATRIX(V,size(r),[rang,1]);
r.L = (S*U')*(r.D*r.L);
r.D = speye(rang);
r.DLmasse={};
r.L = setximasse(r.L,r.DLmasse);
end