function num = getpos(u,num)
%function numl = getpos(u,numg)
% numg : numerotation globale de cette liste
% numl : numerotation locale d'une liste d'element
if nargin==1
num=1:u.nbelem;
else
[rep,num]=ismember(num,u.numelem);
end