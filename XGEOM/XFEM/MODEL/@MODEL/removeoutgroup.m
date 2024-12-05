function M = removeoutgroup(M)

repelim = getnumgroupelemwithfield(M,'lstype','out');

M.groupelem(repelim)=[];
M.nbgroupelem = length(M.groupelem);