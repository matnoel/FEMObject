function S=setlsenrich(S,k,g)
% function S=setlsenrich(S,k,g)
% k : entier (0 : pas d'enrichissement)
% g : groupes à changer (tous par defaut)
if nargin==2
    g=1:getnbgroupelem(S);
end

for i=1:length(g)
    p=g(i);
elem = setlsenrich(getgroupelem(S,p),k);
S = setgroupelem(S,p,elem);
end


