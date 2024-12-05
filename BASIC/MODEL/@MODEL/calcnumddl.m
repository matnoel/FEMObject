function M = calcnumddl(M)

for i=1:M.nbgroupelem
    M.groupelem{i}=calcnumddl(M.groupelem{i},M.node);
end
