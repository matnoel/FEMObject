function M = removenode(M,numnode,choix)
% function M = removenode(M,numnode,choix)
% supprime des noeuds du MODEL M
% si choix = 'global' (par defaut), numnode : numeros des noeuds a eliminer
% si choix = 'local', numnode : position dans la table des noeuds a eliminer

if nargin==2
    choix='global';
end
switch choix
    case 'global'
        [a,b]=ismember(M.number,numnode);
        numnode=find(a);
    case 'local'
        
    otherwise
        error('mauvais choix')
end

if any(numnode)
    garde = setdiff(1:M.nbnode,numnode);
    M.number = M.number(garde(:));
    M.POINT = M.POINT(garde(:));
    M.nbnode = numel(M.number);
end

