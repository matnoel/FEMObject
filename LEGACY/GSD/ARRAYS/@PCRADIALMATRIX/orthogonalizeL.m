function r = orthogonalizeL(r,tol)

%u = normalizeL(u) ; 
%u = uniqueL(u);
%D = speye(u.m);
%for i=1:u.m
% L{i} = u.L(i);
% for j=1:i-1
% D(i,j) = - prodscal(L{i},L{j}) ; 
% L{i} = L{i}+D(i,j)*L{j};    
% end
%end
if nargin==1
    tol=100*eps;
end
D = expectmtimes(r.L,r.L');

s = eig(full(D));

rang = sum(sqrt(abs(s/max(s)))>tol);

%keyboard
opts.disp=0;

%[U,S,V] = svds(D,rang);
if rang<r.m
    [V,S] = eigs(D,rang,'LM',opts);
else
    [V,S] = eig(full(D));    
end
r.m = rang;

%D = diag(diag(S).^(-1/2))*V';
D = V';
r.L = D * r.L ; 
%DV = diag(diag(S).^(1/2))*V';
DV = V';
r.D = speye(r.m);

if iscell(r)
r.V = multimtimes(DV,r.V);
else
    r.V = MULTIMATRIX(double(r.V)*DV',size(r),[rang,1]);
end
r.DLmasse={};
r.L = setximasse(r.L,r.DLmasse);


