function M=setridges(M,F)
M.ridges=cell(1,getnbgroupelem(F));
for i=1:getnbgroupelem(F)
    Fi = keepgroupelem(F,i);
    Fi = removenodewithoutelem(Fi);
    Fi = keepeleminnode(Fi,M.node);
    M.ridges{i} = Fi;
end

M = removeemptyfaces(M);