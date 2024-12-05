function M = removeoutgroup(M)

repelim = [];
for j=1:M.nbgroupelem
    if strcmp(getlstype(M.groupelem{j}),'out')
        repelim = [repelim,j];
    end
end

M.groupelem(repelim)=[];
M.nbgroupelem = length(M.groupelem);