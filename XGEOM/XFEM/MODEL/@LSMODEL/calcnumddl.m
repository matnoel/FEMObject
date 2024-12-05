function M = calcnumddl(M)
% function M = calcnumddl(M)

M = separateelemwithnbddl(M);
for i=1:getnbgroupelem(M)
    elem = calcnumddlenrich(getgroupelem(M,i),getnode(M),M.ls);
    M = setgroupelem(M,i,elem);
end
