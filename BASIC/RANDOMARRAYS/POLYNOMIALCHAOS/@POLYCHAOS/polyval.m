function D = polyval(PC,x,liste)
% function D = polyval(PC,x,liste)
% evaluation des polynomes du chaos multidim en x
% M variables aleatoires
% x : double N-by-M (N : nombre de realisations)
% ou x : 1-by-M cell de N-by-1 double
%
% D matrice de taille N*(P+1)
% liste : numero des polynomes a evaluer (facultatif (par defaut liste=0:P+1))
%   D de taille n*length(liste)

M=PC.M;

indices=PC.indices;
if nargin==3
    indices=indices(liste+1,:) ;
end

d = polyval(PC.RANDPOLYS,indices,x);

D=d{1};
for k=2:M
    D=D.*d{k};
end

