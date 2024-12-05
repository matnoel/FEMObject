function [M,numelem] = keepeleminnode(M,numnode)
% function [M,numelem] = keepeleminnode(M,numnode)
% garde les elements du MODEL M dont tous les noeuds appartiennent a la liste de noeuds numnode
% numelem : numeros des elements conserves
%
% See also MODEL/findeleminnode, MODEL/keepelem, MODEL/keepelemwithnode

if nargin==1
    numnode = getnumber(M.node);
end

numelem = findeleminnode(M,numnode);
M = keepelem(M,numelem);
M = removeemptygroup(M);

M = applyfunctiontofaces(M,@keepeleminnode,numnode);
