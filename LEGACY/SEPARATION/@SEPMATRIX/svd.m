function A = svd(A,tol)
% function A = svd(A,tol)

if getdim(A)>2
    error('utiliser multisvd')
end

A = gathervectors(A);

W1 = A.F{1};
W2 = A.F{2};
sza=size(A.alpha,2);
D = spdiags(A.alpha',0,sza,sza);

[W1,R1] = qr(full(W1));
[W2,R2] = qr(full(W2));
D = R1*D*R2';

[U,S,V] = svd(D);
s = diag(S);
err=sqrt(1-cumsum(s.^2)/sum(s.^2));
m = find(err<tol);
if isempty(m)
    m=getm(A);
else
    m=min(m);
end
Sm = S(1:m,1:m);
Vm = V(:,1:m);
Um=U(:,1:m);

A.F{1} = W1*Um;
A.F{2} = W2*Vm;
A.alpha = diag(Sm);
A.alpha = A.alpha(:)';
A = splitvectors(A);


