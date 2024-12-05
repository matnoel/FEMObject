function M = keepnode(M,numnode)
% function M = keepnode(M,numnode)
% garde des noeuds du MODEL M
% numnode : numeros des noeuds a conserver

[rep,loc]=ismember(numnode,getnumber(M.node));
M.node = getnode(M.node,loc(rep),'local');
M.nbnode = getnbnode(M.node);

M = applyfunctiontofaces(M,@keepnode,numnode);