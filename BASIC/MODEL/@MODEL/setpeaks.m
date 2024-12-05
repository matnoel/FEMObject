function M=setpeaks(M,F)
M.peaks=cell(1,getnbgroupelem(F));

for i=1:getnbgroupelem(F)
    Fi = keepgroupelem(F,i);
    Fi = removenodewithoutelem(Fi);
    Fi = keepeleminnode(Fi,M.node);
    M.peaks{i} = Fi;
end

M = removeemptyfaces(M);
