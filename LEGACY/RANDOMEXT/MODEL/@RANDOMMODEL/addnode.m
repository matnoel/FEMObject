function M = addnode(M,varargin)
% Ajouter des noeuds a un modele
% function M = addnode(M,node)
% node = [1,x1,y1,z1;2,x2,y2,z2;...]  ou un objet NODE ou un POINT
% les noeuds sont ajoutes aux precedents et donc numerotes M.nbnode+[1,2 ...]
%

if nargin==2
    if isa(varargin{1},'MODEL')
M.node=addnode(M.node,varargin{1}.node);
    else
M.node=addnode(M.node,varargin{:});
    end
M.nbnode = get(M.node,'nbnode');
else
    help MODEL/addnode
end
