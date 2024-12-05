function [numelem,numnode]=findelemwithnode(elem,numnode)
if isa(numnode,'NODE')
   numnode = getnumber(numnode); 
end
[a,b]=ismember(elem.connec,numnode);
numelem=elem.numelem(find(sum(a,2)));
numnode=numnode(nonzeros(unique(b)));