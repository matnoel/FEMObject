function a=mean(pcr)
L = pcr.D * pcr.L;
L=mean(L);
if issparse(pcr.V) | issparse(L)
    pcr.V = sparse(pcr.V);
    L = sparse(L);
end

a=sparsemtimes(pcr.V{1},L(1));
for k=2:size(L,1)
a=a+sparsemtimes(pcr.V{k},L(k));
end    

a=reshape(a,size(pcr));

