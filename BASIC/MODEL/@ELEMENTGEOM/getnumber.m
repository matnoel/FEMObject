function num = getnumber(u,num)
% function numg = getnumber(u,numl)
% numl : numerotation locale d'une liste d'element
% numg : numerotation globale de cette liste

if nargin == 1
    num = u.numelem;
else
    num = u.numelem(num);
end