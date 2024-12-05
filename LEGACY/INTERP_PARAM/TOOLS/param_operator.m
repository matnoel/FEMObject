function OPhiA=param_operator(PhiA)

s=size(PhiA,2);
OPhiA=cell(size(PhiA,1),1);
for k=1:size(PhiA,1)
    OPhiA{k}=spdiags(PhiA(k,:)',0,s,s);
end