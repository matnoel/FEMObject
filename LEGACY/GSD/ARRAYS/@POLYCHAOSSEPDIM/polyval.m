function Hc = polyval(PC,x,indices)
% function Hc = polyval(PC,x,indices)
% evaluation des polynomes du chaos multidim en x
% indices : multi-indice [i1,i2,...,iM] du polynome : h_i1(x_1)*h_i2(x_2) ...
% on peut evaluer n polynomes simultanement
% [i11,i21,...;...;i1n,i2n,...]
%
% x : point a evaluer : x=[x1,x2,...,xM]
% on peut evaluer N points simultanement
% x=[x11,...,xM1;x1N,...,xMN]
%
% Hc : matrice de taille N*n
%      Hc = [h_i11(x11)*...*h_iM1(xM1) , ... , h_i1n(x11)*...*h_iMn(xM1) ;
%              ...
%            h_i11(x1N)*...*h_iMn(xMN), ... ,  h_i1n(x1N)*...*h_iMn(xMN)]

M=PC.M;
hc = polyval(PC.RANDPOLYS,indices,x);

Hc=hc{1};
for k=2:M
    Hc=Hc.*hc{k};
end

