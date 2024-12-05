function nbgroupelem = getnbgroupelemwithparam(M,field,value)

num = getnumgroupelemwithparam(M,field,value);
nbgroupelem = length(M.groupelem(num));
