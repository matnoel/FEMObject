function M = removeelem(M,numelem)
% function M = removeelem(M,numelem)
% supprime des elements du MODEL M
% numelem : numeros des elements a eliminer

M.nbelem=0;
for j=1:M.nbgroupelem
    M.groupelem{j}=removeelem(M.groupelem{j},numelem,'global');
    M.nbelem = M.nbelem + getnbelem(M.groupelem{j});
end

M = applyfunctiontofaces(M,@removeelem,numelem);