function M = concatgroupelem(M)
% concatenation des groupes d'elements identiques
% -----------------------------------------------

for i=1:M.nbgroupelem
    mati = getmaterial(M.groupelem{i});
for j=1:i-1
    matj = getmaterial(M.groupelem{j});
    if isempty(mati) & isempty(matj)
        okmat=1;
    else
        okmat = (mati == matj);
    end
    
    okls = lsdatacmp(M.groupelem{j},M.groupelem{i});
    
if strcmp(class(M.groupelem{i}),class(M.groupelem{j})) & okmat & okls
    M.groupelem{j}=concat(M.groupelem{j},M.groupelem{i});
    M.nbgroupelem = M.nbgroupelem-1;
    break
    end
end
end

M.groupelem = M.groupelem(1:M.nbgroupelem);
M=changeelemnumber(M);

