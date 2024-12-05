function isin = lsisin(elem,ls,node)

if nargin==3
    lsxnode = getconnecvalue(elem,ls,node);
else
    lsxnode=ls;
end

%isin = sum(sign(lsxnode)==-1 | sign(lsxnode)==0,2) == getnbnode(elem);
isin = sum(lsxnode<=0,2)==getnbnode(elem);