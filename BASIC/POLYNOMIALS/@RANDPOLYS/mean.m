function Hm = mean(H,indices)

if nargin==2 && isa(indices,'POLYCHAOS')
indices = getindices(indices) ;
elseif nargin==1 && isa(H,'POLYCHAOS')
indices = getindices(H);
end

Hm=sparse(size(indices,1),H.M);
for k=1:H.M
    Hm(:,k) = mean(H.h{k},indices(:,k));
end
Hm = prod(Hm,2);
