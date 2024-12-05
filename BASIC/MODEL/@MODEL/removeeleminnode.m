function [M,numelem] = removeeleminnode(M,numnode)
% function [M,numelem] = removeeleminnode(M,numnode)
% supprime les elements du MODEL M dont tous les noeuds appartiennent a la liste de noeuds numnode
% numelem : numeros des elements supprimes
%
% See also MODEL/findeleminnode, MODEL/removeelem, MODEL/removeelemwithnode

if nargin==1
    numnode = getnumber(M.node);
end

numelem = findeleminnode(M,numnode);
M = removeelem(M,numelem);
M = removeemptygroup(M);

M = applyfunctiontofaces(M,@removeeleminnode,numnode);
