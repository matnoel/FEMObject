function S = setmaterial(S,mat,groupelem)
% function S = setmaterial(S,mat,groupelem)
% associe le materiau mat aux groupes d'elements groupelem

if nargin==2
    groupelem = 1:S.nbgroupelem;
end

for i=groupelem(:)'
    S.groupelem{i} = setmaterial(S.groupelem{i},mat);
end
