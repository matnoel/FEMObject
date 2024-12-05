function n = getgroupelemdim(M)
n=zeros(1,M.nbgroupelem);
for i=1:M.nbgroupelem
    n(i)=getdim(M.groupelem{i});
end
n=unique(n);