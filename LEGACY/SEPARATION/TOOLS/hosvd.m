function [T,err] = hosvd(X,r)
% function T = hosvd(X,r)
% X : N-way tensor
% r : rank of the HOSVD
% T : TSEPMATRIX
% Tensor decompositions and applications - Kolda, Bader, SIAM, 2009

N=ndims(X);
sz=size(X);

if length(r)==1
    r=r*ones(1,N);
elseif length(r)~=N
    error('not implemented')
end

A=cell(1,N);
for n=1:N
    order=[n,1:n-1,n+1:N];
    y=permute(X,order);
    y=reshape(double(y),sz(n),prod(sz([1:n-1,n+1:N])));
    [A{n},temp1,temp2]=svds(y,r(n));
end

T=TSEPMATRIX(N);

T.alpha=ttm(X,A{1}',1);
for n=2:N
    T.alpha=ttm(T.alpha,A{n}',n);
end

for n=1:N
    for j=1:r(n)
        T.F{n}{j}=A{n}(:,j);
    end
end

normX=norm(X(:));

err=abs((normX^2-norm(T.alpha(:))^2))/normX^2;
fprintf('  final error = %d\n',err)
