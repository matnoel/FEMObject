function a=mean(pcr)
L = pcr.L;
for i=1:length(L)
L{i} = mean(L{i});
end

a=sparsemtimes(pcr.V{1},L{1});
for k=2:size(L,1)
a=a+sparsemtimes(pcr.V{k},L{k});
end    

a=reshape(a,size(pcr));

