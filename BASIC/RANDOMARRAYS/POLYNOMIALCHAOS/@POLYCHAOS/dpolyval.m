function Hc = dpolyval(PC,j,x,liste)
% function Hc = polyval(PC,x,liste)
% evaluation des polynomes du chaos multidim en x
% M variables aleatoires
% x : double n-by-M (n : nombre de realisations)
% ou x : 1-by-M cell de n-by-1 double
%
% Hc matrice de taille n*(P+1)
% liste : numero des polynomes a evaluer (facultatif (par defaut liste=0:P+1))
%   Hc de taille n*length(liste)

M=PC.M;

indices=PC.indices;
if nargin==4
    indices=indices(liste+1,:) ;
end

hc = dpolyval(PC.RANDPOLYS,indices,j,x);

Hc=hc{1};
for k=2:M
    Hc=Hc.*hc{k};
end

