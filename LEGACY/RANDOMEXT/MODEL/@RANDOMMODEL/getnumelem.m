function numelem=getnumelem(M)

numelem=[];
for i=1:M.nbgroupelem
 num=getnumber(M.groupelem{i}); 
 numelem=[numelem;num(:)];
end