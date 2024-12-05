function Hc = polyval(PC,x,liste)
% function Hc = polyval(PC,x,liste)
% evaluation des polynomes du chaos multidim en x
% M variables aleatoires
% x : double n-by-M (n : nombre de realisations)
% ou x : 1-by-M cell de n-by-1 double
%
% Hc M cellules de taille (p+1)
% liste : numero des polynomes a evaluer (facultatif (par defaut liste=0:P+1))
%   Hc de taille n*length(liste)

indices = cell(1,getnbdim(PC));
if nargin<=2
    for i=1:getnbdim(PC)
        indices{i} = 0:getn(PC,i)-1;
    end
else
    for i=1:getnbdim(PC)
        indices{i} = liste(:,i);
    end
end

if isa(x,'double')
    x = mat2cell(x,size(x,1),ones(1,size(x,2)));
end

Hc = cell(1,getnbdim(PC));
for i=1:getnbdim(PC)
    Hc{i} = polyval(getpoly(PC,i),indices{i},x{i});
end


