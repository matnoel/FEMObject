function [M,numnode,numelem] = intersect(M,B,varargin)
% function [M,numnode,numelem] = intersect(M,B,varargin)
% function [M,numnode,numelem] = intersect(M,B,'strict',strict)
% intersection d'un MODEL M avec B
% B : MODEL du meme type que M ou GEOMOBJECT
% On cherche l'ensemble numnode des noeuds de M dans B (fonction MODEL/ispointin ou GEOMOBJECT/ispointin)
% puis on garde les elements de M :
% - dont tous les noeuds appartiennent a l'ensemble numnode (fonction MODEL/findeleminnode) si strict est different de 0 (strict ~= 0)
% - dont au moins un noeud appartient a l'ensemble numnode (fonction MODEL/findelemwithnode) si strict est egal a 0 (strict == 0)
% Par defaut, l'option strict est prise egale a 1 pour rester coherent dans le cas ou M et B sont imbriques
% numnode : numeros des noeuds de M appartenant a B
% numelem : numeros des elements de M appartenant a B (au sens strict si strict est different de 0)
%
% See also MODEL/intersect, LSMODEL/setdiff, MODEL/setdiff, MODEL/union

Mini = M;
[M.MODEL,numnode,numelem] = intersect(M.MODEL,B);
if getnblevelsets(M)>0
    M = setlevelsets(M,restrict(getlevelsets(M),Mini,M));
end

if getnblevelsets(M)>0 && isa(B,'MODEL')
    warning('on ne tient compte que des levelset du premier modele')
end

