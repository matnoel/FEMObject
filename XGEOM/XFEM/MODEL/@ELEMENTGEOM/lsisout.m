function isout = lsisout(elem,ls,node)

if nargin==3
    lsxnode = getconnecvalue(elem,ls,node);
else
    lsxnode=ls;
end

%isout = sum(sign(lsxnode)==1 | sign(lsxnode)==0,2) == getnbnode(elem);
isout = sum(lsxnode>=0,2)==getnbnode(elem);