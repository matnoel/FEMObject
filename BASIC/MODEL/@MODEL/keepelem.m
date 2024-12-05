function M = keepelem(M,numelem)
% function M = keepelem(M,numelem)
% garde des elements du MODEL M
% numelem : numeros des elements a conserver

M.nbelem = 0;
for j=1:M.nbgroupelem
    M.groupelem{j}=getelem(M.groupelem{j},numelem,'global');
    M.nbelem = M.nbelem + getnbelem(M.groupelem{j});
end

M = applyfunctiontofaces(M,@keepelem,numelem);