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

indices = cell(1,getnbgroups(PC));
if nargin<=2
    for i=1:getnbgroups(PC)
        indices{i} = 0:(length(PC.PCgroups{i})-1);
    end
else
    warning('prevoir ca proprement')
    for i=1:getnbdim(PC)
        indices{i} = liste(:,i);
    end
end

nd = zeros(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
    nd(i) = length(getgroup(PC,i));
end
if isa(x,'double')
    x = mat2cell(x,size(x,1),nd);
end


Hc = cell(1,getnbgroups(PC));
for i=1:getnbgroups(PC)
    Hc{i} = polyval(PC.PCgroups{i},x{i},indices{i});
end


