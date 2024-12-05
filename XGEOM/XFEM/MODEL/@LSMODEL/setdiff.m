function [M,numnode,numelem] = setdiff(M,B,varargin)
% function [M,numnode,numelem] = setdiff(M,B,varargin)
% function [M,numnode,numelem] = setdiff(M,B,'strict',strict)
% soustraction de B au MODEL M
% B : MODEL du meme type que M ou GEOMOBJECT
% On cherche l'ensemble des noeuds de M dans B (fonction MODEL/ispointin ou GEOMOBJECT/ispointin)
% puis on supprime les elements de M :
% - dont au moins un noeud appartient a cet ensemble (fonction MODEL/removeelemwithnode) si strict est different de 0 (strict ~= 0)
% - dont tous les noeuds appartiennent a cet ensemble (fonction MODEL/removeeleminnode) si strict est egal a 0 (strict == 0)
% setdiff(M,B,'strict',1) correspond au complementaire de intersect(M,B,'strict',0) dans M
% setdiff(M,B,'strict',0) correspond au complementaire de intersect(M,B,'strict',1) dans M
% Par defaut, l'option strict est prise egale a 0 pour rester coherent dans le cas ou M et B sont imbriques
% numnode : numeros des noeuds de M n'appartenant pas a B
% numelem : numeros des elements de M n'appartenant pas a B (au sens strict si strict est different de 0)
%
% See also MODEL/setdiff, LSMODEL/intersect, MODEL/intersect, MODEL/union

Mini = M;
[M.MODEL,numnode,numelem] = setdiff(M.MODEL,B);
if getnblevelsets(M)>0
    M = setlevelsets(M,restrict(getlevelsets(M),Mini,M));
end

if getnblevelsets(M)>0 && isa(B,'MODEL')
    warning('on ne tient compte que des levelset du premier modele')
end

