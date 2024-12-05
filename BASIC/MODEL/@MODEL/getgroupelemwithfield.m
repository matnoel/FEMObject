function groupelem = getgroupelemwithfield(M,field,value)

num = getnumgroupelemwithfield(M,field,value);
groupelem = M.groupelem(num);
