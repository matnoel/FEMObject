function [numelem,numnode]=findeleminnode(elem,numnode)
if isa(numnode,'NODE')
   numnode = getnumber(numnode); 
end
[a,b]=ismember(elem.connec,numnode);
nume = find(sum(a,2)==size(elem.connec,2));
numelem=elem.numelem(nume);
numnode=numnode(nonzeros(unique(b(nume,:))));

