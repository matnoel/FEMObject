function [numelem,numnode] = findeleminnode(M,node)
if isa(node,'NODE')
    node=getnumber(node);
end
numelem = [];
numnode = [];
        for j=1:M.nbgroupelem
          [numelemtemp,numntemp] = findeleminnode(M.groupelem{j},node);
          numelem=[numelem;numelemtemp];
          numnode=union(numnode,numntemp);
        end
