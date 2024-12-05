function bpc = project(apc,pc2,repm)
% function bpc = project(apc,pc2,repm)
% projete a PCARRAY apc on the POLYCHAOS pc2
% if apc is of dimension M1 and pc2 of dimension M2
% repm must be a vector of indices indexing the dimensions of M1
% corresponding to M2

if nargin==2
   repm=1:getM(pc2); 
end

pc1 = apc.POLYCHAOS;

ind2=zeros(getP(pc2)+1,getM(pc1)+1);
ind2(:,[repm,getM(pc1)+1]) = getindices(pc2) ;
[temp,ia]=ismember(ind2,getindices(pc1),'rows');

n=size(apc);n=n(1:end-1);
bpc=myzeros(prod(n),getP(pc2)+1);
apctemp=reshape(apc.MYDOUBLE,prod(n),getP(pc1)+1);
bpc=apctemp(:,ia);
bpc=reshape(bpc,[n,getP(pc2)+1]);
bpc=PCARRAY(bpc,pc2);
