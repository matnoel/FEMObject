function S = separateelemwithnbddlfaces(S)
% function S =  separateelemwithnbddlfaces(S)
% fait en sorte que chaque tous les elements d'un groupe 
% d'elements aient le meme nombre de ddl


for i=1:getnbfacets(S)
S = setfacet(S,i,separateelemwithnbddl(getfacet(S,i)));
end
for i=1:getnbridges(S)
S = setridge(S,i,separateelemwithnbddl(getridge(S,i))); 
end
