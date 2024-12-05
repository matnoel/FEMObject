function [M,numelem] = removeelemwithnode(M,numnode)
% function [M,numelem] = removeelemwithnode(M,numnode)
% supprime les elements du MODEL M dont au moins un noeud appartient a la liste de noeuds numnode
% numelem : numeros des elements supprimes
%
% See also MODEL/findelemwithnode, MODEL/removeelem, MODEL/removeeleminnode

if nargin==1
    numnode = getnumber(M.node);
end

numelem = findelemwithnode(M,numnode);
M = removeelem(M,numelem);
M = removeemptygroup(M);

M = applyfunctiontofaces(M,@removeelemwithnode,numnode);
