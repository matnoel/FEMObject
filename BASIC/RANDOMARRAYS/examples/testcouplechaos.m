p=3;
m=3;
pc = POLYCHAOS(2+m,[p,p,ones(1,m),p+1]);
pc1 = POLYCHAOS(2,p,1)
pc2 = POLYCHAOS(m,1,1)

ind=pc.indices;
ind1=pc1.indices;
ind2=pc2.indices;

vpc = [1:size(ind1,1)];
dpc = [10:10:10*(m+1)];


[t1,loc1]=ismember(ind(:,1:2),ind1(:,1:2),'rows');
[t2,loc2]=ismember(ind(:,3:end-1),ind2(:,1:end-1),'rows');
rep=find(t1 & t2);
apc=zeros(1,size(ind,1));
apc(rep) = vpc(loc1(rep)).*dpc(loc2(rep));
