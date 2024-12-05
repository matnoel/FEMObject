function P = calc_P_facet(M,numfacet)
% function P = calc_P_facet(M,numfacet)
% M : MODEL
% calcul de la matrice permettant de passer des ddl de M1 aux ddl de la
% face numface de M (groupe d'element de M.faces)
% si v1 est un vecteur de ddl defini sur M1
% v2 = P*v1 est le vecteur de ddl correspondant defini sur M2


facet = getfacet(M,numfacet);

numddlmodel = getnumddl(getnode(M));
numddlfacet = getnumddl(getnode(facet));

[locface,locmodel] = ismember(numddlfacet,numddlmodel);

P = sparse(length(numddlfacet),length(numddlmodel));
P(locface,locmodel) = speye(length(locmodel));



