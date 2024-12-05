function iscutorbicut = lsiscutorbicut(elem,ls1,ls2,node)

if nargin==4
ls1xnode = getconnecvalue(elem,ls1,node);
ls2xnode = getconnecvalue(elem,ls2,node);
else
ls1xnode=ls1;
ls2xnode=ls2;
end

iscut1 = lsiscut(elem,ls1xnode);
iscut2 = lsiscut(elem,ls1xnode);
iscutorbicut = iscut1 & ~lsisout(elem,ls2xnode);
isbicut = iscut1 & iscut2;
repbicut = find(isbicut);

for i=1:length(repbicut)
    e=repbicut(i);
    l1 = ls1xnode(e,:)';
    l2 = ls2xnode(e,:)';
    iscutorbicut(e) = islsintersect(l1,l2);
    
end
    