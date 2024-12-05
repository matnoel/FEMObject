function [M,numnode,numelem] = intersect(M,B,varargin)
% function [M,numnode,numelem] = intersect(M,B,varargin)
% function [M,numnode,numelem] = intersect(M,B,'strict',strict)
% intersection d'un MODEL M avec B
% B : MODEL du meme type que M ou GEOMOBJECT
% On cherche l'ensemble des noeuds de M dans B (fonction MODEL/ispointin ou GEOMOBJECT/ispointin)
% puis on garde les elements de M :
% - dont tous les noeuds appartiennent a cet ensemble (fonction MODEL/keepeleminnode) si strict est different de 0 (strict ~= 0)
% - dont au moins un noeud appartient a cet ensemble (fonction MODEL/keepelemwithnode) si strict est egal a 0 (strict == 0)
% intersect(M,B,'strict',0) correspond au complementaire de setdiff(M,B,'strict',1) dans M
% intersect(M,B,'strict',1) correspond au complementaire de setdiff(M,B,'strict',0) dans M
% Par defaut, l'option strict est prise egale a 1 pour rester coherent dans le cas ou M et B sont imbriques
% numnode : numeros des noeuds de M appartenant a B
% numelem : numeros des elements de M appartenant a B (au sens strict si strict est different de 0)
%
% See also MODEL/setdiff, MODEL/union

strict = getcharin('strict',varargin,1);

if isa(B,'GEOMOBJECT') || isa(B,'MODEL')
    repnode = ispointin(B,POINT(M.node));
    numnode = getnumber(M.node,repnode);
    if strict
        [M,numelem] = keepeleminnode(M,numnode);
    else
        [M,numelem] = keepelemwithnode(M,numnode);
    end

% On cherche les noeuds de M et B ayant la meme position geometrique (fonction POINT/intersect)
% puis on repere les elements de B contenant ces noeuds. On garde ensuite dans M les
% elements coincidant avec ces elements (fonction ELEMENTGEOM/intersect)
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
%                 M.groupelem{j} = intersect(M.groupelem{j},B.groupelem{k});
%             end
%         end
%         M.nbelem=M.nbelem+getnbelem(M.groupelem{j});
%     end
else
    error('B must be a GEOMOBJECT or a MODEL')
end

M = removeemptygroup(M);
numnodeelem=getnumnodeelem(M);
M = keepnode(M,numnodeelem);
M = calc_connec(M);
