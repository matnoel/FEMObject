function bpc = prolonge(apc,pc2,repm)
% prolonge a PCARRAY apc on the POLYCHAOS pc2
% if apc is of dimension M1 and pc2 of dimension M2
% repm must be a vector of indices indexing the dimensions of M2
% corresponding to M1

if nargin==2
    repm=1:getM(pc2); 
end
pc1 = apc.POLYCHAOS;

ind1=zeros(getP(pc1)+1,getM(pc2)+1);
ind1(:,[repm,getM(pc2)+1])=getindices(pc1);
[temp,ia]=ismember(ind1,getindices(pc2),'rows');

n=size(apc);n=n(1:end-1);
bpc=myzeros(prod(n),getP(pc2)+1);

apctemp=reshape(apc.MYDOUBLE,prod(n),getP(pc1)+1);
bpc(:,ia)=apctemp;
bpc=reshape(bpc,[n,getP(pc2)+1]);
bpc=PCARRAY(bpc,pc2);