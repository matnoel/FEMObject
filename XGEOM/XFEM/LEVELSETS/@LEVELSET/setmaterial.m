function ls = setmaterial(ls,mat)
% function ls = setmaterial(ls,mat)

if ~isa(mat,'MATERIAL')
    error('rentrer un material')
end
ls.material = mat;
