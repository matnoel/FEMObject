function [v,j]=normsort(v)
v.L = v.D*v.L;
v.D = speye(size(v.D));
L = double(v.L);
Ln = sqrt(sum(L.^2,2));
s=size(v.V);
V=double(v.V);
Vn = sqrt(sum(V.^2,1));

[i,j]=sort(full(Ln(:).*Vn(:)),1,'descend');
v.L = PCMATRIX(L(j,:),[size(L,1),1],getPC(v));
v.V = MULTIMATRIX(V(:,j),s);
