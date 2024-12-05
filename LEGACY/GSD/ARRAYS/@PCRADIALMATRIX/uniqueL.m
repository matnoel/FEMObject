function r=uniqueL(r,tol)
% function r=uniqueV(r,tol)

if nargin==1
    tol=100*eps;
end

D = expectmtimes(r.L,r.L');

s = eig(full(D));
rang = sum(sqrt(abs(s/max(s)))>tol);

%keyboard
opts.disp=0;
if rang<r.m
r.m = rang;
%[U,S,V] = svds(D,rang);
[V,S] = eigs(D,rang,'LM',opts);
%D = diag(diag(S).^(-1/2))*V';
D = V';
r.L = D * r.L ; 
%DV = diag(diag(S).^(1/2))*V';
DV = V';

r.V = MULTIMATRIX(double(r.V)*DV',size(r),[rang,1]);
r.DLmasse={};
r.L = setximasse(r.L,r.DLmasse);
r.D = speye(r.m);


end

