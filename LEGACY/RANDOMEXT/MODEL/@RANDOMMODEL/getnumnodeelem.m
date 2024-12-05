function numnode=getnumnodeelem(M)

numnode=[];
for i=1:M.nbgroupelem
 num=getconnec(M.groupelem{i}); 
 numnode=union(numnode,unique(num));
end
