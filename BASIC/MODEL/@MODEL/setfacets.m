function M=setfacets(M,F)
M.facets=cell(1,getnbgroupelem(F));

for i=1:getnbgroupelem(F)
    Fi = keepgroupelem(F,i);
    Fi = setoption(Fi,'FACE');
    Fi = removenodewithoutelem(Fi);
    Fi = keepeleminnode(Fi,M.node);
    M.facets{i} = Fi;
end

M = removeemptyfaces(M);

