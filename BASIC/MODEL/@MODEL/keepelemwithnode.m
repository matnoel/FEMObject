function [M,numelem] = keepelemwithnode(M,numnode)
% function [M,numelem] = keepelemwithnode(M,numnode)
% garde les elements du MODEL M dont au moins un noeud appartient a la liste de noeuds numnode
% numelem : numeros des elements conserves
%
% See also MODEL/findelemwithnode, MODEL/keepelem, MODEL/keepeleminnode

if nargin==1
    numnode = getnumber(M.node);
end

numelem = findelemwithnode(M,numnode);
M = keepelem(M,numelem);
M = removeemptygroup(M);

M = applyfunctiontofaces(M,@keepelemwithnode,numnode);
