function nbgroupelem = getnbgroupelemwithfield(M,field,value)

num = getnumgroupelemwithfield(M,field,value);
nbgroupelem = length(M.groupelem(num));
