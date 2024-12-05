function A=freemultimatrix(BC,A)
n=length(BC.ddlfree);
rep=repmat(BC.ddlfree(:),1,n) + sqrt(size(A,1))*ones(n,1)*(BC.ddlfree(:)'-1);
A=A(rep(:),:);
