function Indices = calc_multiindices(M,n)
% function Indices = calc_multiindices(M,n)

if length(n)==1
    n = repmat(n,1,M);
end
nmax=max(n);
Indices = zeros(nmax^M,M);
for k=1:M
    for j=1:nmax
        for l=1:nmax^(k-1)
            Indices(l+(j-1)*(nmax)^(k-1):(nmax)^(k):end,k)=j;
        end
    end
end
for k=1:M
    rep=find(Indices(:,k)>n(k));  
    Indices(rep,:)=[];
end
