function us = randomlimit(u,N,Nmax)
% function us = randomlimit(u,N,Nmax)


us=zeros(numel(u),N);
if nargin<3
    if isa(u,'PCMATRIX')
        Nmax=min(N,ceil(1e6/length(getPC(u))));
    else
        Nmax=min(N,5e3);    
    end
end

for k=1:floor(N/Nmax)
    us(:,(k-1)*Nmax+[1:Nmax])=double(random(u,1,Nmax));
%pourcentage(k,floor(N/Nmax))    
end
us(:,floor(N/Nmax)*Nmax+1:end)=double(random(u,1,N-floor(N/Nmax)*Nmax));

