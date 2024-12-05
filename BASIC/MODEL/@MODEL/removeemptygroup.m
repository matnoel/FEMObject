function M = removeemptygroup(M)
% supprime les groupes d'elements vides du MODEL M

repelim = [];
for j=1:M.nbgroupelem
    if isempty(M.groupelem{j}) || getnbelem(M.groupelem{j})==0
        repelim = [repelim,j];
    end
end

M.groupelem(repelim)=[];
M.nbgroupelem = length(M.groupelem);