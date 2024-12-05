function mat = MATERIALS(S)

mat = cell(0,1);
for k=1:S.nbgroupelem
    mat{k} = getmaterial(S.groupelem{k});
end

mat = MATERIALS(mat);