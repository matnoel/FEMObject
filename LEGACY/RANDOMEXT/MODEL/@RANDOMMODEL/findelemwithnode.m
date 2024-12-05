function [numelem,numnode] = findelemwithnode(M,node)
if isa(node,'NODE')
    node=getnumber(node);
end
numelem = [];
numnode = [];
        for j=1:M.nbgroupelem
          [numelemtemp,numntemp] = findelemwithnode(M.groupelem{j},node);
          numelem=[numelem;numelemtemp];
          numnode=union(numnode,numntemp);
        end
