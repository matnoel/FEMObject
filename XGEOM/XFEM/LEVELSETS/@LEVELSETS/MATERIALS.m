function mat = MATERIALS(ls)
% function mat = MATERIALS(ls)

mat = cell(0,1);
for k=1:ls.n
    mat{k} = getmaterial(ls.LS{k});
end
mat = MATERIALS(mat{:});