function valp = getmatrix(A,rep1,rep2)

if nargin==1
for k=1:length(A)    
A{k}=A{k}(:);
end
valp = [A{:}];

else
n1 = length(rep1);
n2 = length(rep2);
valp = sparse(n1*n2,0);
for k=1:length(A)
temp = A{k}(rep1,rep2);
valp = [valp,temp(:)];
end
end

return
