function isbicut = lsisbicut(elem,ls1,ls2,node)
% function isbicut = lsisbicut(elem,ls1,ls2,node)
% il faudrait l'appeler lsisintersect

if nargin==4
ls1xnode = getconnecvalue(elem,ls1,node);
ls2xnode = getconnecvalue(elem,ls2,node);
else
ls1xnode=ls1;
ls2xnode=ls2;
end

isbicut = lsiscut(elem,ls1xnode) & lsiscut(elem,ls2xnode);
repbicut = find(isbicut);

for i=1:length(repbicut)
    e=repbicut(i);
    l1 = ls1xnode(e,:)';
    l2 = ls2xnode(e,:)';
    
    isbicut(e) = islsintersect(l1,l2);
 
end
    