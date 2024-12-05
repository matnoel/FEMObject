function M = removenodewithoutelem(M)
% function M = removenodewithoutelem(M)
% supprime les noeuds isoles (i.e. n'appartenant a aucun groupe d'elements) du MODEL M

co=[];
for i=1:M.nbgroupelem
    co = union(co,unique(getconnec(M.groupelem{i})));    
end
M = keepnode(M,co);

