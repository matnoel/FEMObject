function groupelem = getgroupelemwithparam(M,field,value)

num = getnumgroupelemwithparam(M,field,value);
groupelem = M.groupelem(num);
