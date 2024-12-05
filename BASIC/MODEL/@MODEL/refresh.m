function M= refresh(M)
% M= refresh(M,option)
% finalise le MODEL M
% -> renumerotation des noeuds de 1 a M.nbnode
% -> calcul des connectivites des noeuds et elements
% -> eliminination des noeuds sans element

M = unique(M);
M = calc_connec(M);
% ELIMINATION DES NOEUDS SANS ELEMENT
repnodeelim = find(sum(M.connec.node2elem,1)==0);
if any(repnodeelim)
    M = removenode(M,getnumber(M.node,repnodeelim));
end
M = changenodenumber(M,[1:M.nbnode]);
M = changeelemnumber(M);
M = calc_connec(M);
