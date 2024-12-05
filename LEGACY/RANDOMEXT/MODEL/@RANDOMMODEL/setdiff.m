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
% See also MODEL/setdiff, RANDOMMODEL/intersect, MODEL/intersect, RANDOMMODEL/union

strict = getcharin('strict',varargin,0);

if isa(B,'GEOMOBJECT') || isa(B,'MODEL')
    repnode = ispointin(B,POINT(M.node));
    numnode = getnumber(M.node,repnode);
    if strict
        [M,numelem] = removeelemwithnode(M,numnode);
    else
        [M,numelem] = removeeleminnode(M,numnode);
    end
    
% On cherche les noeuds de M et B ayant la meme position geometrique (fonction POINT/intersect)
% puis on repere les elements de B contenant ces noeuds. On supprime ensuite dans M les
% elements coincidant avec ces elements (fonction ELEMENTGEOM/setdiff)
% elseif isa(B,'MODEL')
%     [P,rep1,rep2] = intersect(POINT(M.node),POINT(B.node));
%     
%     numnode = getnumber(M.node,rep1);
%     numnodeB = getnumber(B.node,rep2);
%     numelem = findeleminnode(M,numnode);
%     numelemB = findeleminnode(B,numnodeB);
%     B = keepelem(B,numelemB);
%     B = keepnode(B,getnumnodeelem(B));
%     B = changenodenumber(B,numnodeB,numnode);
%     
%     M.nbelem = 0;
%     numelem = [];
%     for j=1:M.nbgroupelem
%         for k=1:B.nbgroupelem
%             if isa(M.groupelem{j},class(B.groupelem{k}));
%                 M.groupelem{j} = setdiff(M.groupelem{j},B.groupelem{k});
%             end
%         end
%         M.nbelem=M.nbelem+getnbelem(M.groupelem{j});
%     end
else
    error('B must be a GEOMOBJECT or a MODEL')
end

M = removeemptygroup(M);
numnode=setdiff(getnumber(M.node),numnode);
numnodeelem=getnumnodeelem(M);
M = keepnode(M,numnodeelem);
M = calc_connec(M);
if ~isempty(M.ls)
    M.ls = restrict(M.ls,M);
end

